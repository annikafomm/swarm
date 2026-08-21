# SWARM — Developer Guide

This document describes the *technical* architecture of SWARM: how the frontend and backend fit
together, how users are identified without login, and how uploaded datasets are processed,
registered, and recovered after a restart. For the scientific data format (what goes into
`anndata.obs`/`.obsm`/`.uns` for each score type), see [README.md](README.md).

## Stack overview

```
Angular 19 (frontend/)  <-->  FastAPI + uvicorn (backend/main.py)  -->  Python/R analysis scripts
        nginx (prod only, reverse proxy)                                  (scanpy, squidpy, Seurat,
                                                                            GENIE3, SPONGE, Tangram, ...)
```

- **Frontend**: Angular 19 app in [frontend/](frontend/). Dev server via `npm start` (port 4200,
  proxies `/api` to `localhost:3000` via [frontend/proxy.conf.json](frontend/proxy.conf.json)).
  Production build is served by nginx ([nginx.conf](nginx.conf)), which proxies `/api/*` and a
  fixed list of legacy non-prefixed routes (`/create_session`, `/read_adata`, ...) to the backend
  container, and serves everything else as the Angular SPA.
- **Backend**: FastAPI app in [backend/main.py](backend/main.py) (~2700 lines — the single API
  surface). Delegates heavy lifting to [backend/app.py](backend/app.py)
  (`calculate_scores_helper`, turns a JSON payload into calls to the Python/R pipeline scripts)
  and to [backend/dataset_management.py](backend/dataset_management.py) /
  [backend/dataset_structure.py](backend/dataset_structure.py) (the dataset registry and typed
  `Dataset`/`Params` model classes).
- **Analysis scripts**: `backend/calc_python_scores/`, `backend/calc_R_scores/`,
  `backend/calc_multiome_scores/`, `backend/xenium/`, `backend/dgea/` etc. are invoked as
  subprocesses (`subprocess.run(["python3", ...])` / `["Rscript", ...]`) — not imported — so they
  run in whichever conda env has the right Python/R packages (see `docker_environment.yml`,
  `tangram.yml`).

## Running it locally

- **Devcontainer** (`.devcontainer/`): reinstalls the full conda/R environment
  ([Dockerfile.backend](Dockerfile.backend)) — this takes a long time and blows away installed
  VS Code extensions. Only rebuild it when the environment spec itself changed
  (`docker_environment.yml`, `tangram.yml`, R packages). For everyday code changes, just restart
  the backend/frontend processes inside the existing container.
- **Backend**: `python backend/main.py` (reload is controlled by `UVICORN_RELOAD`, set in
  `.devcontainer/devcontainer.json`). Needs `backend/.env` with a `SESSION_SECRET_KEY` — copy
  `backend/.env.example` and fill in a generated secret; it's gitignored and never baked into the
  Docker image (`.dockerignore`), so every deployment (devcontainer, prod) has its own.
- **Frontend**: `cd frontend && npm start`, or `bash launch.sh` to start both at once.
- **Docker Compose** ([docker-compose.yml](docker-compose.yml)): production-style setup. Backend
  data volumes: `backend/data` is a read-only bind mount (builtin datasets), `swarm-uploads` is a
  named Docker volume holding user uploads + the dataset registry, so it survives container
  recreation.

## User identity — no accounts, cookie-scoped sessions

There is no login. Every "user" is really just an anonymous, server-issued session:

1. On app start, the frontend calls `POST /create_session` ([session.service.ts](frontend/src/app/session.service.ts)).
2. The backend ([main.py](backend/main.py) `create_session`) checks for an existing valid session
   cookie. If found, it's reused as-is (returning user keeps seeing their own uploads). Otherwise
   it mints a fresh `uuid4()`, stores a `SessionData(username=str(uuid))` in an in-memory
   `InMemoryBackend[UUID, SessionData]`, and attaches an httponly, `secure`, `samesite=lax` cookie
   (`SessionCookie`, cookie name `"cookie"|`) via `fastapi_sessions`.
3. The response only ever exposes `username` (= the UUID string) — never the raw session id. The
   frontend cannot choose its own username; that used to be possible and was a dataset-ownership
   bypass (`owner == session_data.username` is the *only* check gating access to uploaded data).
4. Every protected endpoint depends on `Depends(cookie)` + `session_data: SessionData =
   Depends(verifier)` (`BasicVerifier` — sessions never expire server-side beyond in-memory
   backend lifetime; there's a `cleanup_expired_sessions` background task and a manual
   `POST /reset_session` / `POST /delete_session`).
5. Ownership check: an uploaded dataset's `user` field must equal `session_data.username`, or be
   the literal string `"__shared__"` (visible to everyone). Builtin datasets (registered from
   `backend/data/`) always have `user="builtin"` and are always visible
   (`_get_visible_dataset_paths` / `get_datasets` in main.py).

Because sessions live only in memory, **restarting the backend logs everyone out** and forgets
which uploads belong to whom in RAM — but the files and their `job_*_config.json` metadata are
still on disk. That gap is what the "unregistered datasets" recovery flow (below) exists for.

## Data upload flow

1. Frontend posts `multipart/form-data` to `POST /api/upload` with the spatial `.h5ad` plus
   optional single-cell/network/multiome files and ~60 `Form(...)` options (normalization,
   Tangram, network inference, squidpy scores, chromVAR, LIANA+, ...).
2. Backend ([main.py](backend/main.py) `upload()`):
   - **Creates the job directory immediately**, as the very first thing the handler does —
     `backend/uploads/job_<epoch_ms>_<sanitized_username>/` ([main.py:1078-1080](backend/main.py#L1078-L1080)),
     before a single file is saved or validated. If something fails later in the request (a bad
     file, a pipeline crash), the directory is still left behind on disk; see the note on
     `job_*_config.json` below for what does/doesn't get written in that case.
   - Saves every uploaded file into it (`save_file`).
   - If multiome + an `.rds` file was given, converts it to `.h5ad` via
     `Rscript rds_to_h5ad.R` first.
   - Builds a typed `UploadRequest` (Pydantic, [backend/models.py](backend/models.py)) from the
     flat form fields, then a `Params` object ([dataset_structure.py](backend/dataset_structure.py))
     for the OOP dataset model.
   - Calls `calculate_scores_helper(job_dir, payload)` ([backend/app.py](backend/app.py)), which:
     - Translates the JSON payload into CLI flags (`dict2params`) for three independent
       pipelines: a Python script (`calc_python_scores/calc_scores.py`), an R script
       (`calc_R_scores/calc_scores.R`), and optionally multiome-specific Python/R scripts.
     - Runs them sequentially via `subprocess.run(..., check=True)`, writing a `calc_scores.log`
       into a randomly-named scratch subdirectory.
     - For Xenium datasets without Tangram, additionally runs `xenium/compute_backmapping.py` to
       map grid-level scores back onto original cells.
   - Picks the "best" output `.h5ad` in priority order (Xenium cells-with-scores → ST scores →
     Tangram scores → original input) and generates a GeoJSON (`visium_to_geojson.py`) for the
     hexagon/cell map.
   - Wraps everything into a `Dataset` object (`DatasetFactory.create_dataset`, concrete type
     `VisiumDataset`/`XeniumDataset`/`MultiomeDataset`) tagged `user=<session username>`, and
     registers it in the in-memory `DatasetRegistry` (`register_uploaded_dataset`).
   - **Persists a `job_<id>_config.json`** in the job directory containing the full payload —
     this file is the source of truth for recovery, independent of the in-memory registry. It's
     written as the last step of the handler, but since `calculate_scores_helper` swallows its own
     subprocess errors (returns `None` instead of raising), it still gets written even when the
     analysis pipeline fails partway through. Only an exception *earlier* than that (a file-save
     error, bad form data) skips it — leaving a directory with raw input files and no config,
     which the recovery flow below can't do anything with.
3. `DatasetRegistry` ([backend/dataset_management.py](backend/dataset_management.py)) is now purely
   in-memory — `_load_registry`/`_save_registry` and their `dataset_registry.json` disk round-trip
   are commented out (not deleted), along with every call site and the `to_registry_dict`/
   `from_registry_dict` (de)serialization helpers in
   [dataset_structure.py](backend/dataset_structure.py) that only existed to support it. This
   wasn't a behavior change: `_load_registry()` was always called with `clear_uploads=True` (nothing
   in the codebase ever passed `False`), so the `"uploaded"` half of `dataset_registry.json` was
   being written on every mutation but never read back — pure write overhead for a file nothing
   consumed. Uploaded datasets have always been session-scoped and re-registered after a restart via
   the recovery flow below, never repopulated from that JSON file.

## Builtin dataset loading

Builtins were never read from `dataset_registry.json` in the first place — they come from a fixed
list, `DEFAULT_DATASETS` in [main.py](backend/main.py#L198-L211), walked once in the FastAPI
`lifespan` startup hook ([main.py:223-275](backend/main.py#L223-L275)) on **every** process start
(disk persistence is disabled entirely now, see above, but this was already true before that):

- If `<dir>/config.json` exists, `DatasetRegistry.register_builtin_from_config()`
  ([dataset_management.py:643](backend/dataset_management.py#L643)) reads it — note this is
  literally an upload's `job_*_config.json` shape, copied over; e.g. `backend/data/visual_cortex/`
  still has the original job's config sitting alongside the canonical `config.json` used here.
  Paths are resolved via `_resolve_config_path()` (absolute-if-exists → relative to the config's
  own directory → `rglob` by filename as a last resort), which is why builtins tolerate being
  moved around and uploaded datasets (no such fallback) don't.
- Otherwise (`backend/data/` itself has no `config.json`) it falls back to inferring files
  directly by globbing (`adata.h5ad` / `*scores*.h5ad`, `*genie*.csv`, `*sponge*.csv`).
- Either way, the resulting `Dataset` is registered with `user="builtin"` (always visible to
  everyone) — the builtin registration call is what populates the in-memory `datasets["builtin"]`
  dict for this run; nothing round-trips through JSON for it.
- `geojson_path` stored on a builtin dataset is just the URL string `/api/geojson/<id>`, not a
  file path — `GET /api/geojson/{dataset_id}` ([main.py:1896](backend/main.py#L1896)) resolves the
  actual `hexagons.geojson` file per-request via its own fallback chain.

Registry-JSON round-tripping (both builtin and uploaded) is now fully disabled, but it's worth
knowing why, since the commented-out code for it is still in the file for reference: it used to go
through `DatasetFactory.from_registry_dict()`, which is actually broken
(`Dataset.from_registry_dict` → `DatasetFactory.from_registry_dict` → back to
`Dataset.from_registry_dict`, infinite recursion — no subclass ever overrode the base method the
factory delegates to). It was harmless for builtins only because the startup hook immediately
overwrote the failed load with a fresh one. Fixing the recursion (making the base method actually
reconstruct the object instead of delegating — field names line up 1:1 between `to_dict()` and each
concrete class's `__init__`, so a single non-delegating implementation on `Dataset` would work for
all three subclasses) was considered and rejected: even fixed, that path would be strictly worse
than the recovery flow below — it never persisted `params` (the original `Params` object isn't in
`to_dict()` at all) and never re-validates that the referenced output files still exist on disk,
both of which the config-file-based recovery flow does. So the whole mechanism was commented out
rather than repaired.

## Data "reload" / recovery flow

Because the registry drops uploaded datasets on restart (and sessions are in-memory only), a
dataset that was fully processed before a restart would otherwise become invisible even though
its files are untouched on disk. The recovery path:

1. On every app load, once the session is ready, the frontend calls
   `GET /api/unregistered_datasets` ([app.component.ts](frontend/src/app/app.component.ts)).
2. Backend (`DatasetRegistry.get_unregistered_datasets`):
   - Scans `backend/uploads/**/job_*_config.json` (`find_config_files`) — every job directory that
     has a config file is a candidate, registered or not (re-registration is allowed, to support
     recovery after the in-memory copy was cleared). This directory scan is the *only* mechanism
     for bringing uploaded datasets back after a restart — there used to also be a silent
     `_refresh_uploaded_datasets()` step here that replayed `dataset_registry.json` straight back
     into memory, but it went through the same broken `DatasetFactory.from_registry_dict()`
     recursion described above and would silently drop anything it tried to reload. It was removed
     outright rather than fixed: for the small, public user base this app expects, one extra click
     in the recovery dialog below is an acceptable trade for not carrying a second, redundant
     reconstruction path.
   - For each config, infers the type-correct expected output file
     (`_determine_adata_path` — different logic per Visium/Xenium/Multiome and per
     Tangram-enabled/disabled) and reports whether it actually exists on disk
     (`status: "ready" | "missing_files"`).
3. If the list is non-empty, the frontend opens
   `UnregisteredDatasetsDialogComponent`, once per distinct `username` per browser session
   (tracked via `sessionStorage['lastUnregisteredDialogSessionId']` so it doesn't nag on every
   route change).
4. From that dialog the user can:
   - **Register** → `POST /api/register_dataset?dataset_id=...` →
     `DatasetRegistry.register_dataset_from_config` re-reads the job's config JSON, rebuilds the
     concrete `Dataset` object, and registers it under `user=session_data.username` (i.e. the
     *current* session claims it — there is no cross-user ownership recorded in the config file
     itself beyond the directory-name suffix).
   - **Delete** → `POST /api/delete_unregistered_dataset?dataset_id=...` removes the job directory
     entirely.
5. Registering/deleting triggers `DatasetService.loadAvailableDatasets()` again so
   `GET /api/datasets` (builtin + datasets owned by `session_data.username` + `__shared__`) reflects
   the change immediately.

**Practical implication**: uploads are never truly lost on backend restart, but they *are*
tied to whichever session re-registers them — there's no durable, cross-session "my datasets"
list beyond what's encoded in `backend/uploads/`.

## Which AnnData is "the" AnnData

A processed dataset carries several `.h5ad` files and it is easy to assume one of them is the raw
upload. None of them is — **the raw upload is never retained as a dataset field**. It lives only in
the transient `UploadRequest.files.spatial_h5ad` and is not represented on the `Dataset` classes at
all.

The pipeline calls the same scoring routine twice on two different objects
([calc_scores.py:154-157](backend/calc_python_scores/calc_scores.py#L154-L157)):

- **`adata_st_scores.h5ad`** — the uploaded spatial object (after optional filtering/normalization)
  with all scores computed **directly on the spots**.
- **`adata_tg_scores.h5ad`** — the *Tangram output* object (single-cell reference mapped onto the
  spatial coordinates) with the same scores computed on **that** instead. Only exists if Tangram
  ran. Additionally carries `obsm["tangram_ct_pred"]`, which LIANA+ needs and which chromVAR scores
  are projected through.

**`adata_path` ("primary AnnData") is not a file of its own — it is a pointer** resolved by
`_determine_adata_path()` ([dataset_management.py:439](backend/dataset_management.py#L439)):

| Dataset | `adata_path` resolves to |
|---|---|
| Visium, no Tangram | `adata_st_scores.h5ad` |
| Visium + Tangram | `adata_tg_scores.h5ad` |
| Multiome | always `adata_tg_scores.h5ad` |
| Xenium, no Tangram | `xenium_cells_with_grid_scores.h5ad` |
| Xenium + Tangram | `adata_tg_scores.h5ad` |
| fallback | whatever `output_files["adata_path"]` says |

So for a Multiome dataset, `adata_path`, `adata_tg_scores_path` **and** `tangram_adata_path` all
point at the same file. That is not a bug — `tangram_adata_path` is documented as an *alias* for
`adata_tg_scores_path` ("ignored if `adata_tg_scores_path` is provided",
[dataset_structure.py:822](backend/dataset_structure.py#L822)), i.e. it refers to projected
**output**, not to the single-cell reference input. Do not label it as a "reference" anywhere in the
UI; that mistake has been made and reverted once already.

## How `cell_type` and `leiden` are decided, and what consumes them

These two are routinely conflated but they are produced by completely different mechanisms, and the
asymmetry explains several "why is this tab empty" questions.

### `leiden` — machine-decided, and only as a side effect

Preferred path: the upload already has an `obs` column literally named `leiden`.

Fallback: SWARM computes it, but only when **three** conditions hold simultaneously
([calc_scores.py:91-96](backend/calc_python_scores/calc_scores.py#L91-L96)) — centrality,
co-occurrence or neighborhood enrichment was requested; **and** `leiden` is not already in `obs`;
**and** that score's cluster-key argument is still the literal string `"leiden"`. Only then does
`clustering()` ([preprocessing_functions.py:41](backend/preprocessing/preprocessing_functions.py#L41))
run `sc.pp.neighbors` → `sc.tl.umap` → `sc.tl.leiden(flavor="igraph", n_iterations=2,
directed=False)`, at scanpy's default resolution. The in-code comment is explicit that this is a
*"makeshift solution for when no cluster key is provided"* and deliberately not user-configurable.

Two consequences:

- **Select no squidpy cluster score and no clustering ever happens.** `builtin_heart_multiome` is
  exactly this case (`scores.squidpy = false` in its config): there is no `leiden` column at all,
  so `meta.leiden_cluster_annotations` comes out `{}` and Cluster Information / Co-occurrence /
  Neighborhood Enrichment are legitimately unavailable. That is correct behavior, not breakage.
- **The per-score cluster-key parameters are a trap.** The upload form exposes `cluster_cs`,
  `cluster_co` and `cluster_nhood` separately, and squidpy honours them — writing
  `uns["<yourcolumn>_centrality_scores"]` etc. But
  [visium_to_geojson.py](backend/visium_to_geojson.py) only ever reads the three hardcoded
  `leiden_*` names, so anything computed under a different key is silently unreachable from the UI.
  Labels must additionally cast to `int` (`.astype(int)`), or the whole annotation block empties.

### `cell_type` — biology you have to bring

Supplied path: an `obs` column. Tangram's cluster mode expects it under the name `cell_type`
(`--cell_label`, default `cell_type`).

Derived path (Tangram), two distinct stages that are easy to conflate:

1. `map_cells_to_space(..., cluster_label=cell_label)` produces the mapping matrix, then
   `project_cell_annotations` writes **`obsm["tangram_ct_pred"]`** — a *composition*: one
   continuous score per cell type per spot. Not a label.
2. `add_dissociated_annotations_to_obs()`
   ([calc_tangram.py:347](backend/calc_python_scores/calc_tangram.py#L347)) takes
   `ct_pred.idxmax(axis=1)` → **`obs["cell_type_dissociated"]`**, the argmax winner, plus one
   `prob_<celltype>_dissociated` column per type and a `uns["dissociated_prob_column_map"]` index.

Mapping mode is `'clusters'` by default and forced to `'cells'` for multiome
([calc_tangram.py:240-243](backend/calc_python_scores/calc_tangram.py#L240-L243)).

Two non-obvious behaviors:

- **Tangram never overwrites `cell_type`.** It writes `cell_type_dissociated` alongside it, so a
  Tangram run leaves you with three related things: the original annotation, the composition
  matrix, and the derived argmax label.
- **Tangram post-processing *requires* a pre-existing `cell_type` column.** `mask_obs_key` defaults
  to `'cell_type'` and `add_dissociated_annotations_to_obs` raises `KeyError` if it is absent, so
  Tangram cannot bootstrap cell types from nothing — it needs a seed annotation on the spatial side.

#### Can Tangram run without any cell-type annotation?

Not as currently wired. Five stages touch `cell_label`; three are avoidable and two are not:

| Stage | Needs an annotation? | Avoidable? |
|---|---|---|
| `select_genes` | `ctg` and `spapros` do, `hvg` does not; selecting no mode returns `None` so Tangram uses all overlapping genes | yes |
| `map_cells_to_space(cluster_label=…)` | only in `mode='clusters'`, which aggregates cells *by* the label | yes — multiome forces `mode='cells'`, where Tangram ignores `cluster_label` |
| `project_genes` | mirrors the mode above | yes, same way |
| `project_cell_annotations(annotation=cell_label)` ([calc_tangram.py:283](backend/calc_python_scores/calc_tangram.py#L283)) | yes — called unconditionally, outside any mode check | **no** |
| `add_dissociated_annotations_to_obs(mask_obs_key="cell_type")` ([calc_tangram.py:469](backend/calc_python_scores/calc_tangram.py#L469)) | yes — raises `KeyError` when the column is missing | **no** |

Two wrinkles in those last two rows:

- They want the annotation on **different objects**. `project_cell_annotations` reads `cell_label`
  from `ad_map.obs` (the single-cell reference), while `mask_obs_key` is checked against
  `ad_ge.obs` (the spatial side). A Tangram run therefore needs a cell-type column on *both* inputs.
- `mask_obs_key` is hardcoded to the string literal `"cell_type"` rather than `args.cell_label`, so
  even `--cell_label annotation_final` still demands a column named exactly `cell_type`. This looks
  unintended.

Via the web app it is stricter again: **`cell_label` has no form field**, so it is always
`cell_type` and cannot be pointed at a differently-named column without invoking the script directly.

Conceptually an annotation-free multiome run *should* be possible: the mapping matrix
`adata_map` (FoPra's M) is what chromVAR and footprint projection consume, and M does not need cell
types in `cells` mode — so spatial chromVAR would still work. Only the cell-type composition
products would be lost (`tangram_ct_pred`, `cell_type_dissociated`, and hence LIANA+
cell-composition × TF activity). Supporting it would mean guarding `project_cell_annotations` on
the label being present, making `mask_obs_key` follow `cell_label`, and allowing
`add_dissociated_annotations_to_obs` to be skipped.

### What consumes which

| Consumer | `leiden` | `cell_type` |
|---|---|---|
| Squidpy cluster stats → `meta.leiden_cluster_annotations` | yes | — |
| Sidebar tabs | Cluster Information, Co-occurrence | Cell Information |
| Map palette | `leidenColorScale` (Tableau10) | `colorScale` (Set2) |
| Cluster outlining (`extendCluster`) | leiden view only | — |
| LIANA+ cell-composition × TF activity | — | `liana.composition_column`, default `tangram_ct_pred` |
| Differential motif activity | — | grouping column (`spot_groupby`, hardcoded) |
| Footprints | — | grouping (`cluster_by`, defaults to `cell_type`) |
| DGEA | selectable grouping | **default** grouping (`selectedDgeaObsCol`) |

`cell_type` is also what the map falls back *to*: the initial view priority is `regulatory_scores`
→ `cell_type` → first available property, and `getAvailableView()` hardcodes `'cell_type'` as its
fallback argument.

The two categorical palettes are deliberately different so it is obvious at a glance whether the
map is showing clusters or cell types — they used to share one ordinal scale.

**Summary of the asymmetry**: `leiden` is structure *discovered from* the data — cheap and
automatic, but it only materialises if a squidpy score happened to ask for it, and it is only
readable under one exact name. `cell_type` is biology *brought to* the data — nothing derives it
from scratch, several downstream analyses hardcode it as their grouping, and it is the map's
fallback view. This is why a dataset with no clustering at all (heart) still renders sensibly: it
has `cell_type`, and the ATAC tabs group by `cell_type` rather than `leiden`.

## Frontend Spatial Map, Cell, Cluster, and Tab Synchronization

The spatial visualization ([HexagonPlotComponent](frontend/src/app/hexagon-plot/hexagon-plot.component.ts) and [HexagonViewComponent](frontend/src/app/hexagon-view/hexagon-view.component.ts)) follows strict synchronization rules between the map canvas, color views, sidebar tabs, and cell/cluster selection state:

### 1. Cell Clicking Behavior
- **When colored by "Leiden" (`selectedView === 'leiden'`)**:
  - **Cluster selection**: Selects the clicked cell's entire Leiden cluster (`selectedCluster = clusterId`).
  - **Map borders**: Outlines all cells belonging to that Leiden cluster with a thin black border (`1.4px`), and outlines the clicked cell itself with a thick black border (`3px`).
  - **Dimming**: Non-cluster cells are dimmed to `0.6` opacity; cluster member cells are rendered at `1.0` opacity.
  - **Sidebar tab**: Jumps directly to the **"Cluster Information"** tab.
  - **Background updates**: Updates the Cell Information panel in the background so switching to Cell Info displays the clicked cell's metadata immediately without losing context.
- **When colored by ANY other view (e.g. `cell_type`, `gene_expression`, `regulatory_scores`, `co_occurrence`, etc.)**:
  - **Single-cell selection**: Selects only that specific cell (`selectedCell = cell`).
  - **Map borders**: Outlines *only* the clicked hexagon with a thick black border (`3px`).
  - **No cluster borders or dimming**: Cluster outlines and background dimming are strictly disabled (`opacity: 0.8` for all unselected cells).
  - **Sidebar tab**: Jumps directly to the **"Cell Information"** tab.
  - **Background updates**: Cluster statistics and neighborhood enrichment heatmaps are calculated in the background without switching tabs or applying cluster outlines.

### 2. "Color hexagons by" Dropdown Changes
- Changing the dropdown property (`onColorbyPropertyChange`) automatically switches the sidebar to the corresponding tab:
  - `leiden` $\rightarrow$ **Cluster Information**
  - `cell_type` / numeric metadata $\rightarrow$ **Cell Information**
  - `gene_expression` $\rightarrow$ **Gene Expression**
  - `regulatory_scores` (or specific Genie3/SPONGE scores) $\rightarrow$ **Regulatory Scores**
  - `co_occurrence` $\rightarrow$ **Cluster Information** (Spatial Co-occurrence Matrix)
  - `ligand_receptor_*` / LIANA+ $\rightarrow$ **Ligand-Receptor Relationships**
  - `cell_comp_tf_activity_*` $\rightarrow$ **Cell Composition TF Activity**
  - `tf_activity` $\rightarrow$ **TF Activity**
  - `pathway_activity` $\rightarrow$ **Pathway Activity**
  - `dgea` $\rightarrow$ **DGEA**
  - `chromvar_*` / spatial correlation $\rightarrow$ **ChromVar spatial correlation : Moran's I / Geary's C**
  - `diff_motif_*` $\rightarrow$ **Differential Motif Activity**
  - `footprints` $\rightarrow$ **Footprints**
  - `grn_evaluation` $\rightarrow$ **GRN Evaluation**
  - Any other dataset attribute $\rightarrow$ **Further Attributes**
- **Border & Dimming Cleanup**: If the new color view is *not* `leiden`, cluster extension outlines and dimming are cleared immediately (`activeClusterId = null`, standard `0.8` opacity across the canvas).

### 3. Cluster Information Tab & Dropdown
- **Switching to "Cluster Information" Tab**:
  - Automatically switches the map color view to `leiden` (`selectedView = 'leiden'`).
  - Highlights and outlines the active cluster (or auto-selects default cluster 0 if no cluster is active yet).
  - Keeps `selectedCell = null` if entering cold (displays empty state guide in Cell Info).
- **Selecting a Cluster in the Cluster Info Dropdown (`selectCluster`)**:
  - Switches map color view to `leiden`.
  - Outlines the newly selected cluster on the map.
  - Resets single-cell thick border (`selectedCell = null`) if the previously selected cell was not part of the new cluster.

### 4. Further Attributes Tab & Spatial Metrics Engine
- **Dynamically Discovered Dataset Attributes**:
  - Automatically scans any dataset (Visium, Xenium, or custom GeoJSON) for features/properties not mapped to dedicated tabs.
  - **Unique Cell Identifiers**: Attributes where each cell has a unique ID ($N_{\text{unique}} \approx N$, e.g. `barcode`, `observation_joinid`, `centroid`) are excluded from the visual table and presented in the *Cell Information* tab.
  - **Dataset-wide Constant Metadata**: Properties with 1 unique value across all spots ($N_{\text{unique}} = 1$, e.g. `Sample_ID`, `treatment`, `assay`, `organism`, `disease`, `tissue`) are categorized under the dedicated **Sample Information** category in *Cell Information*.
  - **Visualizable Attributes**: Cell composition fractions, region annotations, density scores, and custom numerical/categorical columns appear in the *Further Attributes* table.
- **On-The-Fly Spatial Statistics**:
  - Computes **Min**, **Max**, **Mean**, **Variance ($s^2$)**, **Moran's I ($I$)**, and **Getis-Ord General $G$ ($G$ & $Z$-score)** in real-time from the spatial $k=6$ neighbor graph.
  - Selecting any attribute recolors the hexagons on the spatial map with appropriate continuous or categorical scales.

### 5. Cell Properties & Missing Values
- **Dataset-wide empty properties**: Properties that have no non-empty values across any cell in the entire dataset are filtered out completely.
- **Cell-level missing values**: If a property exists in the dataset but has a null, undefined, empty string, or NaN value for the selected cell, the property key is displayed with a styled `NaN` badge.

### 5. Canvas Aspect Ratio & Legend Bounding
- **Square Projection**: Canvas uses a `1000x1000` SVG viewBox fitted with `d3.geoIdentity().fitExtent([[20, 20], [980, 980]], featureCollection)` inside a 1:1 aspect ratio container (`#hexbin`), providing equal margin spacing on all four sides.
- **Floating Legend**: Positioned at `(x: 20, y: 20)` within SVG bounds with positive internal offsets and drop shadow, preventing any clipping at canvas edges.
