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
2. **Spatial Autocorrelation**: Projection smooths expression across neighboring spots, naturally yielding higher Moran's $I$ values than raw spatial measurements.
3. **`cell_type` vs. `leiden`**:
   * `cell_type`: Required for Tangram deconvolution, differential motif grouping, and footprinting.
   * `leiden`: Computed automatically by Scanpy only when spatial clustering metrics are requested and no pre-existing cluster column exists.

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
2. **Leiden View (`selectedView === 'leiden'`)**:
   * Clicking a spot selects its parent cluster (`activeClusterId = clusterId`).
   * Cluster member spots receive a **thin black border (`1.4px`)**.
   * The clicked spot receives a **thick black border (`3px`)**.
   * Non-member spots have a **transparent border (`1px`)**.
   * Navigates directly to the **Cluster Information** tab.
3. **Feature / Continuous Views (`cell_type`, `gene_expression`, `regulatory_scores`, etc.)**:
   * Clicking a spot outlines *only* that individual hexagon with a **thick black border (`3px`)**.
   * All other spots have a **transparent border (`1px`)** (no cluster outlines).
   * Navigates directly to the **Cell Information** tab.
4. **Color-By Dropdown**:
   * Changing the active feature property synchronizes the sidebar to the corresponding tab and clears active cluster outlines (`activeClusterId = null`).
5. **Compare View**:
   * Renders main and compare datasets in synchronized half-height containers with independent color-scales and feature selectors.

---

## 9. Known Gotchas & Troubleshooting

* **Dataset `config.json` Network Paths**: Ensure `output_files` contains `genie_network_path` and `sponge_network_path` pointing to the filtered network CSVs; otherwise, interactive D3 subgraphs (`/geneset_connections_*`) will return empty results.
* **Named AnnData `obs` Index in LIANA+**: When unpacking loadings, `liana_bivariate.py` uses `_index_by_first_column()` to prevent `KeyError` on datasets where `adata.obs.index.name` is not `"index"`.
* **Case-Sensitivity in Mouse Datasets**: LIANA+ consensus resources use uppercase human symbols. When processing mouse datasets, run backfills with `--keep-var-case` to prevent permanent mutation of `adata.var_names`.
* **Corrupted `uns['neighbors']`**: Datasets processed without Scanpy clustering may contain partial neighbor dictionaries. `backfill_squidpy_scores.py` automatically clears malformed `uns['neighbors']` prior to computing graphs.
