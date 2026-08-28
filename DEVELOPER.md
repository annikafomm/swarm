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
