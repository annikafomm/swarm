![all_yellow_cursive_white](all_yellow_cursive_white.svg)

# SWARM

This is the Github-Repository of the MoPiTas Project in the Master Practical course of DaisyBioLab 2025.

## Dataset-Table

https://docs.google.com/spreadsheets/d/1RH_lB4OEUEgZMIYu59_RTwcPmkC-ciyApEnlhvPSJ7E/edit?usp=sharing 

## Overview

![Overview](overview.png)

## Before Starting SWARM
### Example Data
Use the download_example_data.sh Script to autimatically download and move the example data to the right places.

## Backend

Start the backend development server:

```sh
python backend/main.py
```

## Frontend

Start the frontend development server:

```sh
cd frontend
npm start
```
OR use
```sh
bash launch.sh
```
## Format Requirements
### Visium

The web-tool is working with the anndata-format. This is a tutorial on how counts and scores should be saved.

- **Counts**
  - Stored in `anndata.X`
- **Observations (per Spot)**
  - Stored in `anndata.obs` (all observation columns will be shown, so make sure to clean your anndata before upload). 
- **Clustering**
  - Save the Leiden clustering for Leiden-based scores in `anndata.obs["leiden"]`
- **Gene-wise scores**
  - Moran's I: stored in `anndata.uns["moranI"]` as DataFrame with columns:
    - `I`
    - `pval_norm`
    - `var_norm`
    - `pval_norm_fdr_bh`
  - Geary's C: stored in `anndata.uns["gearyC"]` as DataFrame
    - `C`
    - `pval_norm`
    - `var_norm`
    - `pval_norm_fdr_bh`
- **Leiden cluster scores**
  - Leiden centrality: stored in `anndata.uns["leiden_centrality_scores"]` as DataFrame with columns: 
    - `degree_centrality`
    - `average_clustering`
    - `closeness_centrality`
  - Leiden co-occurrence: stored in `anndata.uns["leiden_co_occurrence"]` as dictionary with entries:
    - `intervals`: interval boundaries, array of shape `(n_intervals - 1,)`
    - `occ`: 3D array of shape `(n_clusters, n_clusters, n_intervals)` (Make sure to select the right num_intervals in Upload form)
  - Leiden neighborhood enrichment: stored in `anndata.uns["leiden_nhood_enrichment"]` as dictionary with entries:
    - `counts`: array of shape `(n_clust, n_clust)`
    - `zscore`: array of shape `(n_clust, nclust)`
- **Regulon scores**
  - Gene sets: stored in `anndata.uns["genie_genesets"]` and `anndata.uns["sponge_genesets"]` as dictionaries with format:
    - `regulator`: `[gene1, gene2, ..., geneN]`
  - AUCell: stored in `anndata.obsm["aucell_scores_{sponge|genie3}"]` as DataFrames with regulators as columns.
  - GSVA: stored in `anndata.obsm["spongeffects_GSVA_scores_{sponge|genie3}"]` as DataFrames with regulators as columns.
  - ssGSEA: stored in `anndata.obsm["spongeffects_ssGSEA_scores_{sponge|genie3}"]` as DataFrames with regulators as columns.
  - Viper: stored in `anndata.obsm["viper_scores_genie3"]` as DataFrame with regulators as columns.
- **DGEA**
  DGEA is currently implemented as a Seurat-based workflow.  
The input AnnData object is first exported into a Seurat-compatible format, then differential expression is computed in R using `Seurat::FindMarkers`, and the resulting DGEA results are merged into the dataset metadata for frontend visualization.

### Workflow

1. Export the AnnData object to Seurat-compatible files:
   - `matrix.mtx`
   - `genes.tsv`
   - `barcodes.tsv`
   - `metadata.csv`

2. Load the exported data in Seurat and create a Seurat object.

3. Optionally filter the dataset to `in_tissue == 1` if the column is available.

4. Normalize the Seurat object.

5. Run differential expression analysis across selected metadata columns.

6. Write the results to `dgea_results.json`.

7. Merge the DGEA results into the GeoJSON metadata so they can be accessed by the frontend.

### Groupings

The current DGEA workflow is configured to run on the following metadata columns if they are available:

- `cell_type`
- `leiden`
- `niche_cluster`

For each grouping column, the workflow computes:

- all pairwise group comparisons
- one-vs-all comparisons
- an aggregate marker view based on the best adjusted p-value observed across comparisons

### Statistical settings

DGEA is computed with `Seurat::FindMarkers` using the Wilcoxon test.

Current default settings include:

- `min_cells = 3`
- `top_n = 200`
- `min_pct = 0.10`
- `logfc_threshold = 0.25`

For the DGEA heatmap context, the workflow additionally uses:

- `heatmap_n_top = 10`
- `heatmap_padj_cutoff = 0.05`
- `heatmap_logfc_cutoff = 0.5`
- normalized expression values from the Seurat `"data"` slot

### Output structure

The DGEA results are stored as JSON under `meta["dgea"]`.

For each grouping column, the output contains:

- available levels
- group sizes
- aggregated marker results
- individual comparisons

Each comparison stores:

- comparison ID
- group names
- group sizes
- skipped status and skip reason if applicable
- a marker table
- a heatmap context containing scaled and raw expression values across groups

Comparisons are skipped automatically if one of the groups does not satisfy the minimum cell requirement.

### Frontend integration

The frontend reads DGEA data from the dataset metadata and provides a dedicated **DGEA** tab.

Users can:

- select the grouping column
- choose group 1 and group 2
- switch to a one-vs-all comparison
- inspect the marker table for the selected comparison
- view a heatmap of the selected comparison
- project selected genes back onto the main plot via **Show on plot**

- **Tangram**
  - Cell type compositions: stored in `anndata.obsm["tangram_ct_pred"]` as DataFrame with cell types as columns.
- **LIANA+ scores**
  - Ligand receptor relationships:
    - Cosine similarity: stored in `anndata.obsm["ligand_receptor_cosine_similarity"]` as array of shape (n_spots, n_interactions)
    - P-value: stored in `anndata.obsm["ligand_receptor_p_value"]` as array of shape (n_spots, n_interactions)
    - Category: stored in `anndata.obsm["ligand_receptor_category"]` as array of shape (n_spots, n_interactions)
    - NMF factors: stored in `anndata.obsm["ligand_receptor_NMF_factors"]` as DataFrame with factors as columns.
    - Global scores: stored in `anndata.uns["ligand_receptor_global_scores"]` as DataFrame with columns:
      - `cosine_similarity_mean`
      - `cosine_similarity_std`
      - `ligand_receptor_morans`: Moran's R score
    - Interaction names: stored in `anndata.uns["liana_columns"]["ligand_receptor"]` as array of shape (n_interactions,)
  - Cell type composition - TF activity similarity:
    - Cosine similarity: stored in `anndata.obsm["cell_comp_tf_activity_cosine_similarity"]` as array of shape (n_spots, n_interactions)
    - Category: stored in `anndata.obsm["cell_comp_tf_activity_category"]` as array of shape (n_spots, n_interactions)
    - Global scores: stored in `anndata.uns["cell_comp_tf_activity_global_scores"]` as DataFrame with columns:
      - `cosine_similarity_mean`
      - `cosine_similarity_std`
    - Interaction names: stored in `anndata.uns["liana_columns"]["cell_comp_tf_activity"]` as array of shape (n_interactions,)
  - TF activity:
    - ULM score: stored in `adata.obsm["tf_activity_score_ulm"]` as DataFrame with TFs as columns.
    - ULM p-adjusted value: stored in `adata.obsm["tf_activity_padj_ulm"]` as DataFrame with TFs as columns.
  - Pathway activity:
    - MLM score: stored in `adata.obsm["pathway_activity_score_mlm"]` as DataFrame with pathways as columns.
    - MLM p-adjusted value: stored in `adata.obsm["pathway_activity_padj_mlm"]` as DataFrame with pathways as columns.

We furthermore use GRNs computed using sponge and genie3, both CSV-formatted.

Genie:

```csv
regulatoryGene,targetGene,weight
FOXM1,KIF20A,0.06364686682258
```

Sponge:

```csv
geneA,geneB,df,cor,pcor,mscor,p.val,p.adj
ENSG00000182141,ENSG00000258630,1,0.121684631317227,0.0764637550413806,0.0452208762758468,0.029426,0.519790968019617
```

### Xenium-specific

## Overview

For Xenium data, the pipeline does not run all analyses directly on the original cell-level object.  
Instead, it first creates a **grid-level representation** and performs downstream spatial analyses on that representation.

The overall workflow is:

1. Load the original Xenium AnnData object.
2. Optionally filter and normalize the data.
3. Build a **grid-level AnnData** for downstream scoring.
4. Keep a **cell-level copy** of the original Xenium data.
5. Map each original Xenium cell to its nearest grid spot.
6. Run downstream spatial analyses on the **grid / spot level**.
7. Optionally run Tangram on the Xenium grid-level object.
8. If Tangram is **not** used, propagate grid-level results back to the original Xenium cells.
9. If Tangram **is** used, keep the results on the **spot / grid level** for visualization.

## Xenium processing

Because Xenium data is available at **single-cell resolution**, the pipeline converts it into a **grid-based spatial representation** before scoring. This makes downstream spatial analyses more practical while still preserving the original cell-level data for later use.

For Xenium datasets, the pipeline:

- creates a grid-level AnnData object
- stores a mapping from original cells to grid spots
- performs downstream analyses on the **spot / grid level**

Important intermediate files are:

- `xenium_map.h5ad`
- `xenium_grid.h5ad`

## Tangram behavior

Tangram is optional.

When Tangram is enabled, the pipeline uses the **grid-level Xenium object** as the spatial input. Tangram is therefore not run directly on the original Xenium cell-level object.

Tangram supports multiple gene-selection modes:

- `ctg`
- `hvg`
- `spapros`

An optional user-defined gene list can also be provided.

## Spatial score calculation

Depending on the selected options, the Xenium workflow can compute:

- spatial neighbors
- LIANA+
- centrality scores
- co-occurrence
- neighborhood enrichment
- Moran’s I
- Geary’s C

The main score outputs are:

- `adata_st_scores.h5ad`
- `adata_tg_scores.h5ad`

## Back-mapping behavior

### Without Tangram

If Xenium is processed **without Tangram**, the pipeline performs a back-mapping step after scoring.  
In this step, grid-level results are propagated back to the original Xenium cells.

The final output is:

- `xenium_cells_with_grid_scores.h5ad`

### With Tangram

If Xenium is processed **with Tangram**, **no back-mapping to the original Xenium cells is performed**.

Instead:

- Tangram-derived results remain on the **spot / grid level**
- the result is visualized on the **spot level**

This avoids the additional cost of propagating Tangram-derived results back to every original Xenium cell.

## Output behavior

After the pipeline finishes, the backend selects output files in this order:

1. `xenium_cells_with_grid_scores.h5ad`
2. `adata_st_scores.h5ad`
3. `adata_tg_scores.h5ad`
4. otherwise the original spatial input

This means:

- **Xenium without Tangram** is usually visualized from `xenium_cells_with_grid_scores.h5ad`
- **Xenium with Tangram** falls back to `adata_tg_scores.h5ad`
- therefore, **Xenium + Tangram is displayed on spot / grid level**

## Frontend visualization behavior

The frontend contains Xenium-specific visualization logic.

### Comparison mode

Comparison mode is **disabled for Xenium datasets**.

When a Xenium dataset is selected:

- comparison mode cannot be toggled on
- an existing comparison view is automatically disabled
- the compare button is disabled in the UI

### Performance handling

For Xenium visualization, the frontend uses dedicated performance handling:

- Xenium datasets are detected either by metadata (`data_type == "xenium"`) or by large feature counts
- a reduced base layer is rendered for performance
- a detail window is used for more focused interactive inspection
- Xenium rendering is handled differently from standard Visium rendering

This means the frontend is optimized for large Xenium datasets and does not treat Xenium exactly like a standard Visium hexagon view.

## Summary

In this pipeline, Xenium is handled through a **grid-first workflow**:

- Xenium cells are converted into a grid representation
- downstream analyses run on the **spot / grid level**
- Tangram uses the **grid-level Xenium object**
- **without Tangram**, results are mapped back to cells
- **with Tangram**, results remain on the **spot / grid level** for visualization



### Multiom-specific

