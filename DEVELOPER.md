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
