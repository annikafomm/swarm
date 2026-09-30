# SWARM — Developer Guide

This document describes the technical architecture of SWARM: frontend/backend integration, dataset lifecycle, session management, analysis pipelines, and map synchronization.

---

## 1. Architecture & Tech Stack

```
Angular 19 (frontend/)  <──HTTP/WS──>  FastAPI + Uvicorn (backend/main.py)  ──Subprocesses──>  Python / R Pipelines
    │                                              │                                             │
    └─ nginx (reverse proxy in prod)               ├─ In-memory DatasetRegistry & LRU Cache      ├─ scanpy, squidpy, liana
                                                   └─ Session Cookie Auth (HMAC signed)          └─ Seurat, GENIE3, SPONGE, Tangram
```

* **Frontend**: Angular 19 SPA (`frontend/`).
  * Dev server: `npm start` (port 4200, proxies `/api` to `localhost:3000` via `proxy.conf.json`).
  * Production: Static bundle served by Nginx (`nginx.conf`).
* **Backend**: FastAPI app (`backend/main.py`).
  * Orchestrates pipelines (`backend/app.py`), dataset models (`backend/dataset_structure.py`), and registry management (`backend/dataset_management.py`).
  * Exposes cached AnnData queries (`_ADATA_CACHE`) and GeoJSON endpoints.
* **Pipelines**: Python and R scripts invoked as external subprocesses (`subprocess.run`) across conda environments.

---

## 2. Local Development & Setup

### Requirements & Config
* Copy `backend/.env.example` to `backend/.env` and set `SESSION_SECRET_KEY`.
* Set `UVICORN_RELOAD=true` for live backend code reloading.

### Commands
```bash
# Start backend
python backend/main.py

# Start frontend
cd frontend && npm start
when error with file watcher: WATCHPACK_POLLING=true npm start

# Start both
bash launch.sh

# Docker Compose (Production emulation)
docker compose up --build
```

## Data Privacy, Legal Compliance & Tracking

SWARM implements the **Reduced Tracking Case** under German and EU data protection law (GDPR / DSGVO, TDDDG, DDG, BayDSG):

- **No Cookie Banner Required**:
  - **Cookieless Matomo**: Matomo is self-hosted on TUM servers (`exbio.wzw.tum.de`) and operates strictly without cookies (`_paq.push(['disableCookies'])`) and with Do-Not-Track enabled (`_paq.push(['setDoNotTrack', true])`). IP addresses are masked by 3 bytes, visitor profiles/logs are disabled, and data is purged after 180 days.
  - **Strictly Essential Session Cookie**: The HttpOnly FastAPI session cookie (`cookie`) is strictly technically necessary (§ 25(2) No. 2 TDDDG) to maintain workspace isolation and on-demand calculations. It does not track users across websites.
- **Legal Footer & Privacy Policy Dialog**:
  - The application includes a legal footer linking the TUM ExBio Impressum (`/impressum#de`), Legal Notice (`/impressum#en`), and an interactive bilingual Privacy Policy dialog (`PrivacyDialogComponent`).
  - The Privacy Policy modal includes an interactive **Matomo Opt-Out switch** so users can exercise their Art. 21 GDPR right to object at any time.

### Devcontainer
* Devcontainer (`.devcontainer/Dockerfile.backend`) builds the full Conda + R environment.
* Only rebuild the container when dependency manifests (`docker_environment.yml`, `tangram.yml`) change.

---

## 3. Session & Authentication Model

SWARM uses anonymous, cookie-based session management (`fastapi_sessions`):

1. **Session Creation**: `POST /create_session` assigns a signed, `HttpOnly`, `SameSite=Lax` session cookie containing a random UUID (`SessionData(username=uuid)`).
2. **Access Control**:
   * **Builtin Datasets**: Tagged `user="builtin"`; accessible to all users.
   * **Uploaded Datasets**: Tagged `user=<session_uuid>` or `user="__shared__"`; accessible only by the session owner.
3. **Session Persistence**: Sessions live in-memory. Restarting the backend invalidates active session IDs but preserves uploaded files on disk. Datasets can be re-registered via `GET /api/unregistered_datasets`.

---

## 4. Dataset Lifecycle & Resolution

### Dataset Hierarchy & Output Files

`DatasetFactory` instantiates typed dataset models based on modalities:

| Dataset Type | Primary AnnData Pointer (`adata_path`) | Key Associated Files |
| :--- | :--- | :--- |
| **Visium (Standard)** | `adata_st_scores.h5ad` | `hexagons.geojson`, `genie_network_filtered_st.csv` |
| **Visium (Tangram)** | `adata_tg_scores.h5ad` | `adata_map.h5ad`, `hexagons.geojson` |
| **Multiome** | `adata_tg_scores.h5ad` | `spot_obj_chromvar.rds`, `adata_map.h5ad`, ATAC fragments |
| **Xenium** | `xenium_cells_with_grid_scores.h5ad` / `adata_tg_scores.h5ad` | Cell geometries, grid mappings |

### AnnData File Resolution
* **Live API Queries** (`_resolve_adata_path` in `main.py`):
  * Multiome: Checks `adata_st_scores_path` $\to$ `adata_tg_scores_path` $\to$ `adata_path`.
  * Xenium: Checks `tangram_adata_path` $\to$ `xenium_grid_adata_path` $\to$ `adata_path`.
  * Visium: Checks `tangram_adata_path` $\to$ `adata_path`.
* **GeoJSON Generation** (`_determine_adata_path` in `dataset_management.py`):
  * Resolves the primary scored AnnData object to serialize spot properties.

### Builtin vs. Uploaded Registration
* **Builtin Datasets** (`backend/data/<dataset_name>/`):
  * Loaded on startup via `DEFAULT_DATASETS` in `main.py`.
  * Configuration defined in `<dataset_name>/config.json`.
  * Paths in `config.json` resolve relative to the dataset directory.
* **Uploaded Datasets** (`backend/uploads/job_<timestamp>_<user>/`):
  * Pipeline writes `job_<timestamp>_config.json` on completion.
  * Registry scans these JSON files on `GET /api/unregistered_datasets` to allow recovery across backend restarts.

---

## 5. Expression Data Semantics (`st` vs. `tg`)

Tangram-enabled runs score the **projected expression object** (`adata_tg_scores.h5ad`), not the measured spatial spots (`adata_st_scores.h5ad`):

```python
# calc_scores.py execution logic
if args.tangram:
    compute_spatial_scores(adata_tangram, "tg", ...)   # Single-cell projected space
else:
    compute_spatial_scores(adata_work, "st", ...)      # Measured spatial spots
```

### Important Implications
1. **Values**: On Tangram datasets, the *Gene Expression* tab, regulatory scores, and spatial metrics reflect imputed/projected gene densities.
2. **Spatial Autocorrelation — corrected**: the line that used to stand here ("projection smooths expression, naturally yielding higher Moran's I") understates what is actually going on. It is not generic smoothing, it is a specific, measurable artifact that can inflate *or destroy* autocorrelation depending on the estimator — see 5.1 below.
3. **`cell_type` vs. `leiden`**:
   * `cell_type`: Required for Tangram deconvolution, differential motif grouping, and footprinting.
   * `leiden`: Computed automatically by Scanpy only when spatial clustering metrics are requested and no pre-existing cluster column exists.

### 5.1 Why the projected object's spatial scores are not trustworthy

Investigated in depth on the heart AVN multiome (3153 spots); full numbers, code and raw output
in [`heart_projection_control_findings.md`](heart_projection_control_findings.md) and
[`backend/scripts/heart_projection_control.py`](backend/scripts/heart_projection_control.py).
Two independent problems compound, not one.

**(a) Tangram's `project_genes` is a linear operator: `X_space = M.T @ X_sc`**
(`tangram/utils.py:368`, confirmed by reading the installed package — no `.raw` reference exists
anywhere in `tangram`, it uses whatever `adata_sc.X` holds). `M` has shape
(reference cells x spots) with **every row summing to exactly 1.0** (sd 5.6e-08) — each reference
cell's probability distribution over spots. Its per-spot column sums (`m[s] = sum_c M[c,s]`,
"mapping mass") have **CV 0.2194**.

That mass factor is the *entire* per-gene spatial signal in the projected object:

| measurement (heart, 26,868 genes) | value |
|---|---:|
| Moran's I of `m` itself | **+0.8525** |
| median \|corr(gene's projected field, `m`)\| | **0.833** (43% of genes > 0.9) |
| median per-gene Moran's I, raw | +0.591 |
| median per-gene Moran's I, after regressing `m` out | **-0.0013** |

Rank-based regulon estimators (AUCell, VIPER, GSVA, ssGSEA) are invariant to any per-spot
rescaling `X[s,:] -> lambda_s * X[s,:]`, so they implicitly divide `m[s]` back out — and once it
is gone, nothing spatially structured is left. Measured regulon Moran's I on the projected
object: AUCell **-0.0014**, VIPER **+0.0063** (both from the pipeline's own stored scores).
Cross-validated with an independent Python implementation: **-0.0084**.

**Practical consequence**: `n_perms=None` (the pipeline's only autocorrelation null, see
Gotcha "Corrupted `uns[\'neighbors\']`" below and `add_to_adata.py`) tests against spatial
*randomness*, which `m` violates before any biology enters — hence up to **98.6%** of genes
looking "significant" on a Tangram-projected object. Any Moran's I computed on such an object
needs a mapping-aware null (permute the reference profile across cells, re-project through the
same `M` — do **not** use an i.i.d. random profile, that under-generates the null by ~30x) or it
is measuring the mapping, not the tissue.

**(b) SWARM feeds `project_genes` log-scale input; the reference pipeline uses linear counts.**
Verified empirically on heart: both `ad_sc.raw` (`calc_tangram.py:236`) and the spatial `ad_sp.X`
are log1p-scale (non-integer, max ~6.9-10).

*Fitting and projecting are not two steps with independently choosable scales — they are the
same linear operator, applied twice, confirmed by reading the installed package.* The fit step's
loss (`mapping_optimizer.py`) computes

```python
M_probs = softmax(self.M, dim=1)
G_pred  = torch.matmul(M_probs.t(), S)          # identical form to project_genes' M.T @ X_sc
gv_term = cosine_similarity(G_pred, G, dim=0).mean()
```

where `S = adata_sc[:, training_genes].X` and `G = adata_sp[:, training_genes].X`
(`mapping_utils.py:260,269`) — **no transform between `.X` and either matrix**. So whatever scale
`adata_sc.X`/`adata_sp.X` hold at call time is the scale the cosine-similarity loss is fit
against, and it is the same scale `project_genes` later multiplies by `M.T`. There is no
"fit in log space, project in linear space" option in the package as shipped — the input scale
choice is made once, before `pp_adatas`, and both steps inherit it.

Given that, feeding log-scale data is the wrong choice specifically *because* of the second use.
Cosine similarity is scale-tolerant for the fit (it depends only on vector direction, so a
monotone per-gene transform mostly reorders which genes dominate the mean, it does not break the
objective). Linear projection is not: `project_genes` computes `sum_c w_c * log(1+x_c)`, but the
mixture Tangram is meant to estimate is `log(1 + sum_c w_c*x_c)`. Since `log` is concave, Jensen's
inequality means the first is systematically <= the second — projected expression is compressed
relative to the true mixture, worst wherever a spot mixes cells with very different expression
(exactly the boundary regions AVN-style analyses target). Fitting tolerates the wrong scale;
projecting does not.

*Checked against primary sources, not summaries — the official tutorial notebook was fetched and
parsed directly (grepped the raw `.ipynb` JSON, not a paraphrase):*
after `sc.read_h5ad`, the very first thing the [official Tangram tutorial](https://github.com/broadinstitute/Tangram/blob/master/tutorial_tangram_without_squidpy.ipynb)
does is verify the loaded data is raw counts (`np.unique(ad_sc.X.toarray()[0,:])`, immediately
following its own stated rule "if the data are in integer format, that probably means they are in
raw count"), applies `sc.pp.normalize_total` only — **zero occurrences of `log1p` or `.scale(` in
the entire notebook** — and passes the spatial object through with **no preprocessing at all**.
Cell-level commentary: *"mapping works great with raw data... [normalize_total is] light
pre-processing."* This is a genuine deviation from documented Tangram usage, not an equally-valid
alternative reading.

**What the papers say — confirmed, not just the code/tutorial.** The original paper
([Biancalani et al. 2021, *Nat. Methods*](https://www.nature.com/articles/s41592-021-01264-7))
states in its Methods that the single-cell and spatial count matrices are **normalized for
library size only** before mapping — no log transform. This matches the tutorial and the
installed package exactly: `S`/`G` in the fit loss and `X_sc` in `project_genes` are all
library-size-normalized counts, never logarithmized, at every point the paper, the code, and the
tutorial notebook can be checked against each other. There is no ambiguity left on this point.

Given that, a later benchmarking paper
([refinement strategies for Tangram, Bioinformatics/ISMB 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12261478/))
is simply wrong on this detail — it states in its own methods that *"As suggested by Tangram's
developers, both datasets are normalized by the total counts over all genes per cell or spot and
**logarithmized** before running Tangram,"* which is not what the original paper or the official
tutorial does. Treat that specific sentence as a mistaken paraphrase, not a competing convention.
The same paper is still useful independent evidence for scoring the measured object over the
projected one, on grounds unrelated to this point: it documents that Tangram gives
**"inconsistent results over repeated runs"** and that **"genes expressed in fewer cells or spots
tend to have poorer predictions"** (i.e. the projected layer's reliability degrades exactly on
the sparser, more spatially-interesting genes) — both orthogonal to, and consistent with, the
mapping-mass finding above.

**Bottom line: SWARM's heart pipeline feeds log1p data into an interface that both the paper and
the official tutorial specify as library-size-normalized-only.** This is confirmed as a bug
against the documented method, not a defensible alternative reading of ambiguous guidance.

**(c) Scoring the measured object instead is correct, but is not a drop-in fix.** Measured
object regulon Moran's I: AUCell **+0.4393**, VIPER **+0.5221**, GSVA **+0.6280**, ssGSEA
**+0.6194** — all bracketing the +0.578 predicted independently in advance. But:
  * The two objects disagree on gene-identifier namespace (projected: symbol index +
    `ensembl_id` column; measured: Ensembl index + `SYMBOL` column) and must be reconciled
    before scoring — see `backend/scripts/prepare_measured_for_scoring.py`.
  * **Rank-based R scorers are pathologically slow on real (sparse) Visium counts, and a
    per-gene detection filter does not fix it.** The projected object is fully dense (0 genes at
    zero); AUCell/VIPER completed on it in 16s/12s. The measured object is ~3% dense: **15,092 of
    33,538 genes are zero in every spot**, and even after filtering to genes detected in >=10
    spots (13,222 genes, 617/634 GENIE3 regulators and 92.6% of targets retained), a *single* spot
    still has ~11,600 zero-valued genes forming one enormous tie group for AUCell/GSVA/ssGSEA to
    rank through (vs. a largest tie group of 7 on the projected object). That is a per-**spot**
    problem, not a per-**gene** one, so gene filtering barely helps it: scoring took 5h+ once
    unfiltered (killed, no output) and ~6h15m once filtered. Budget accordingly; do not expect a
    gene-count filter alone to make this fast.
  * Scores are per-spot x per-regulon (gene-axis independent), so after scoring the smaller
    filtered object they were transplanted back onto the full symbol-indexed 33,538-gene object
    to preserve full Gene Expression coverage — see the same script's usage in the heart config.

**Where this now lives for heart**: `adata_st_scores_path` in `backend/data/heart/config.json`
points at `plasmidpoop/adata_st_scores_measured.h5ad` (regulatory scores computed on the
*measured* object), while `adata_tg_scores_path` is untouched. Because `_resolve_adata_path`
(live API — gene values, obsm tables) prefers `adata_st_scores_path` while `_determine_adata_path`
(registration + GeoJSON) uses `adata_tg_scores_path`, the map's baked-in properties are still
projection-derived but live-fetched regulatory scores now come from the measured object.
`calc_scores.py` also gained an opt-in `-score_measured_too` flag (default off, existing runs
bit-for-bit unchanged) that scores the measured object into the `st` slot alongside a Tangram
run's `tg` slot, so future re-runs do not need this by-hand transplant.

---

## 6. Regulatory Networks (GENIE3 & SPONGE)

SWARM calculates regulatory activities (AUCell, GSVA, ssGSEA, VIPER) by scoring spatial matrices against user-supplied or precomputed GRNs.

### Network File Contracts
* **GENIE3**: CSV with headers `regulatoryGene, targetGene, weight` (matched by **Gene Symbol**).
* **SPONGE Interaction**: CSV with headers `geneA, geneB, df, cor, pcor, mscor, p.val, p.adj` (matched by **Ensembl ID**).
* **SPONGE Centrality**: CSV with headers `gene, degree, eigenvector, betweenness, page_rank` (matched by **Ensembl ID**).

### Scoring Execution Flow
`calc_R_scores/calc_scores.R` parses networks $\to$ builds regulon modules $\to$ computes scores $\to$ `calc_python_scores/add_to_adata.py` merges outputs:
* `obsm`: `aucell_scores_*`, `spongeffects_GSVA_scores_*`, `spongeffects_ssGSEA_scores_*`, `viper_scores_*`.
* `uns`: `genie_genesets` (symbols), `sponge_genesets` (Ensembl), `*_moranI`, `*_gearyC`.

### Backfilling Existing Datasets
Use the CLI backfill utilities in `backend/scripts/` to calculate missing score layers without re-running preprocessing:

```bash
# 1. Backfill Squidpy spatial metrics (Leiden, Moran's I, Co-occurrence, Centrality)
python3 backend/scripts/backfill_squidpy_scores.py --adata <path_to_adata>

# 2. Backfill LIANA+ signaling (LR interactions, TF activities, pathways)
python3 backend/scripts/backfill_liana_scores.py --adata <path_to_adata>

# 3. Backfill Regulatory Scores (AUCell, GSVA, ssGSEA, VIPER)
python3 backend/scripts/backfill_regulatory_scores.py --adata <path_to_adata> \
    --genie3-network backend/networks/GENIE3/heart_gtex_lv/genie3_heart_lv_network.csv \
    --sponge-network backend/networks/SPONGE/heart_gtex_lv/sponge_heart_lv_interactionNetwork.csv \
    --sponge-analysis backend/networks/SPONGE/heart_gtex_lv/sponge_heart_lv_networkAnalysis.csv

# 4. Regenerate GeoJSON to update the frontend
python3 backend/visium_to_geojson.py --adata <path_to_adata> \
    --outpath <dataset_dir>/hexagons.geojson --data_type visium
```

---

## 7. Performance & Resource Throttling

To prevent worker exhaustion during heavy R/Python computations:

1. **Thread Pinning**:
   ```bash
   OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
   ```
   Set in `Dockerfile.backend` and `docker-compose.yml` to prevent BLAS libraries from oversubscribing host CPU cores.
2. **Pipeline Concurrency Gate**:
   `MAX_CONCURRENT_PIPELINES` (default: `2`) throttles execution in `main.py` via `asyncio.Semaphore` to protect host RAM from simultaneous multi-gigabyte AnnData allocations.
3. **Backend AnnData Caching**:
   `_load_adata_cached` maintains an in-memory LRU cache of recently queried AnnData objects to avoid repetitive disk I/O on spot inspections.

---

## 8. Frontend Component Architecture & Map Sync

### Component Decomposition
* **`HexagonPlotComponent`** (`frontend/src/app/hexagon-plot/`): Main coordinator managing state, dataset switching, and toolbar controls.
* **`HexagonViewComponent`** (`frontend/src/app/hexagon-view/`): Encapsulated D3 canvas rendering hexagonal/spot geometry, color interpolations, and SVG legends.
* **Sidebar Panels**: Extracted dedicated components for clean separation of concerns:
  * `CellInfoPanelComponent` / `ClusterInfoPanelComponent`
  * `RegulatoryScoresPanelComponent` / `RegulatoryTablesPanelComponent`
  * `ChromvarCorrelationPanelComponent`
  * `FurtherAttributesPanelComponent`

### Map Selection & Tab Synchronization Rules
1. **Selection is Purely Border-Based (No Dimming)**:
   * Cell fills are never dimmed when selecting cells or clusters; highlighting is achieved entirely through SVG stroke width and color.
2. **Categorical Views** (any property where `isContinuousScale()` is false — Leiden clusters, `cell_type`, or any other categorical `obs` column; integer values with ≤20 uniques, or non-numeric strings, count as categorical):
   * Clicking a spot selects its parent cluster under the *currently active* property (`activeClusterValue` in `HexagonViewComponent`).
   * Cluster member spots receive a **thin black border (`1.4px`)**.
   * The clicked spot receives a **thick black border (`3px`)**.
   * Non-member spots have a **transparent border (`1px`)**.
   * Navigates directly to the **Cluster Information** tab, which shows a **"Cluster by" dropdown** (`HexagonPlotComponent.categoricalProperties()`) to switch which categorical attribute drives clustering. Options exclude the same "silly" attributes (spatial/grid coordinates, IDs, constant columns) that `FurtherAttributesPanelComponent` hides, via the shared `attribute-filters.ts` (`isSpatialOrIdentifierKey`, `isUninformativeAttribute`) — both surfaces stay in sync from one heuristic instead of two.
   * Leiden-only "heavy" data (centrality averages, co-occurrence, neighborhood enrichment — computed only under `leiden_*` uns keys) is hidden by `ClusterInfoPanelComponent`'s `isLeidenView` input, and the corresponding side effects (`updateCoOccurrenceTable`, `renderFootprintPlots`, `updateSubgraphGenie3`, nhood heatmap) are gated behind `view === 'leiden'` in `HexagonPlotComponent` whenever a different categorical property is active.
   * Switching the active clustering property re-derives the selected cluster (`selectClusterForActiveProperty`): the selected cell's value under the new property if a cell is selected, otherwise the same default-cluster heuristic used on initial load (`autoSelectDefaultCluster` — for Leiden, prefers a cluster with co-occurrence data; otherwise the most populous value).
3. **Continuous Views** (`gene_expression`, `regulatory_scores`, etc.):
   * Clicking a spot outlines *only* that individual hexagon with a **thick black border (`3px`)**.
   * All other spots have a **transparent border (`1px`)** (no cluster outlines).
   * Navigates directly to the **Cell Information** tab.
4. **Color-By Dropdown**:
   * Changing the active feature property synchronizes the sidebar to the corresponding tab and clears active cluster outlines (`activeClusterValue = null`) unless the new property is categorical, per rule 2.
5. **Tab-Jump Suppression** (`HexagonPlotComponent.suppressTabJump` / `withoutTabJump()`): a single choke point that both jump primitives (`jumpToTabByLabel`, `jumpToTab`) check before moving the sidebar. Any control that lives inside a tab and triggers a property change that would otherwise re-derive and jump to a (possibly different) tab should wrap its call in `withoutTabJump(() => ...)`, rather than threading a new boolean parameter through every intermediate method — the pattern `onColorbyPropertyChange`'s own `skipTabJump` parameter had to use. Used today by the "Cluster by" dropdown, so switching the clustering property never jumps away from the Cluster Information tab.
6. **Compare View**:
   * Renders main and compare datasets in synchronized half-height containers with independent color-scales and feature selectors.

---

## 9. Known Gotchas & Troubleshooting

* **Dataset `config.json` Network Paths**: Ensure `output_files` contains `genie_network_path` and `sponge_network_path` pointing to the filtered network CSVs; otherwise, interactive D3 subgraphs (`/geneset_connections_*`) will return empty results.
* **Named AnnData `obs` Index in LIANA+**: When unpacking loadings, `liana_bivariate.py` uses `_index_by_first_column()` to prevent `KeyError` on datasets where `adata.obs.index.name` is not `"index"`.
* **Case-Sensitivity in Mouse Datasets**: LIANA+'s `consensus` resource uses uppercase human symbols, so `run_liana` used to do `adata.var.index = adata.var.index.str.upper()` in place — permanently rewriting mouse symbols (`Xkr4` -> `XKR4`) in the object that then gets written to disk. **Fixed**: `calc_liana.run_liana` now takes `organism=` (default: auto-detected from symbol casing) and selects LIANA's mouse-native `mouseconsensus` resource instead of touching the data. On visual_cortex this also *improves* LR coverage (89.5% of resource genes present vs 81.3% for the uppercase-vs-human match). `backfill_liana_scores.py` takes the same `--organism {auto,human,mouse}`. For other species use `li.rs.get_hcop_orthologs()` + `li.rs.translate_resource()` on a *copy* of the resource.
* **Repairing an already-mangled dataset**: `str.upper()` is not invertible — correct mouse casing follows no single rule (`Rb1cc1`, `4732440D04Rik`, `mt-Nd1`, `AI597479` each need a different one), so re-capitalising cannot recover it. Use `scripts/restore_var_symbols.py --adata <file> --source <uncorrupted file>`, which copies symbols positionally and refuses unless the correspondence is provable (same var count, elementwise upper-case match, no upper-case collisions). `var_names` is not the only casualty: every score stage that ran *after* LIANA wrote its results keyed by the mangled symbols, so the script repairs those in the same pass —
  * `uns['moranI']`, `uns['gearyC']` — one row per gene
  * `uns['{viper,aucell,spongeffects_GSVA,spongeffects_ssGSEA}_scores_genie3_{moranI,gearyC}']` — one row per regulon
  * `uns['genie_genesets']` / `uns['sponge_genesets']` — both the regulon keys and the target-gene lists
  * `obsm['*_scores_{genie3,sponge}']` — regulon column names

  Two guards keep it from over-reaching, and both matter in practice: a table is only rewritten if **every** one of its names resolves against the restored var index, and obsm is restricted to the `_genie3`/`_sponge` suffixes the pipeline uses. Without them, `obsm['pathway_activity_score_mlm']` (whose PROGENy column names include a literal `EGFR`) and `obsm['proportions_class']` (whose cell-class labels `Sst`, `Vip`, `Pvalb`, `Lamp5`, `Sncg`, `Meis2` are all real mouse genes) would be silently relabelled into something that means the wrong thing.

  Deliberately **not** repaired: `uns['ligand_receptor_global_scores']` and `uns['liana_columns']['ligand_receptor']`. Those hold `LIGAND^RECEPTOR` names from the **human** `consensus` resource — not mangled mouse symbols, so there is nothing to map them back to. Mouse-native pairs require re-running the LR stage:
  ```bash
  python3 backend/scripts/backfill_liana_scores.py --adata <file> --organism mouse \
      --skip cell_comp_tf_activity pathway_activity --force
  ```
  After any restore, regenerate the GeoJSON — the frontend reads these names from `meta`, not from the h5ad.

* **Human datasets are affected too, just less visibly.** It is tempting to assume `.upper()` is a no-op on human symbols. It is not: HGNC writes uncharacterised open reading frames with a lowercase `orf` (`C1orf159`), and there are a few hundred of them. The heart builtin had 321 such genes upper-cased to `C1ORF159`.

  Heart also shows why "just restore the backup" is not automatically right. Its symbols were **entirely lowercase** before the LIANA backfill (`al627309.1`, `linc01409`, `samd11`) — itself non-standard — so the upper-casing actually *improved* 26,547 of 26,868 symbols and broke only the 321 `orf` ones. Reverting to the backup would have regressed almost everything. The authoritative source was instead `var['SYMBOL']` in `HCAHeartST11290662_adata_final_annotation.h5ad` (proper HGNC casing, 0 upper-case collisions), reached with the script's `--source-column SYMBOL --allow-unmatched`:
  ```bash
  python3 backend/scripts/restore_var_symbols.py \
      --adata backend/data/heart/plasmidpoop/adata_tg_scores.h5ad \
      --source backend/data/heart/HCAHeartST11290662_adata_final_annotation.h5ad \
      --source-column SYMBOL --allow-unmatched --backup
  ```
  `--allow-unmatched` is needed because the source is a different annotation build: 738 of heart's genes are absent from it (`LINC01409`, `PRXL2B`) and are already correct, so they are left alone. That run corrected 317 symbols plus 4,000 rows each in `uns['moranI']` and `uns['gearyC']` — which were *lowercase*, having been computed before the upper-casing, and so disagreed with `var_names` in the opposite direction. One gene, `C12ORF81`, is absent from the source and remains mis-cased; `MORF4L1`/`MORF4L2`/`MORF4L2-AS1` merely contain the substring and are correctly upper-case.

  **Impact on heart's scores was small**, because the human reference networks are almost entirely upper-case: `consensus` has **zero** mixed-case genes, so the ligand-receptor results were unaffected and need no re-run. `progeny(human)` has 170 mixed-case targets of 17,610 and `collectri(human)` 46 of 6,675, so pathway and TF activity missed roughly 1% of their networks while the symbols were upper-cased. Re-running those stages now recovers that 1%.
* **Corrupted `uns['neighbors']`**: Datasets processed without Scanpy clustering may contain partial neighbor dictionaries. `backfill_squidpy_scores.py` automatically clears malformed `uns['neighbors']` prior to computing graphs.
