# Multiome Pipeline – Output Files

All outputs are written inside the job output directory:

```
uploads/job_<id>/<outdir>/
```

where `<outdir>` is a randomly chosen name (`plasmidpoop`, `junkDNA420`, `kackhaufen1`, or `dumpase1`).

---

## Directory overview

```
uploads/job_<id>/
├── <uploaded_multiome.rds>                   ← input (user upload)
├── <fragments.tsv.gz>                        ← input (user upload, optional)
├── <fragments.tsv.gz.tbi>                    ← input (user upload, optional)
└── <outdir>/
    ├── calc_scores.log
    ├── adata_st_scores.h5ad
    ├── adata_tg_scores.h5ad                  ← always written by R multiome pipeline, bc always tangram
    ├── adata_map.h5ad                        ← always written by R multiome pipeline, bc always tangram
    ├── adata_map.X.csv                       ← always written by R multiome pipeline, bc always tangram
    ├── adata_map.var.csv                     ← always written by R multiome pipeline, bc always tangram
    ├── motif_to_tf.csv                       ← if any motif analysis flag is set
    ├── global_motif_analysis.rds             ← always written by R multiome pipeline
    ├── spot_obj_chromvar.rds                 ← if --chromvar + --differential_motif_activity
    ├── spot_obj_footprints.rds               ← if --footprinting
    ├── dissociated_obj_footprints.rds        ← if --footprinting
    ├── footprint_<motif_id>.pdf MIGHT CHANGE TO     footprint_<cell_type>_<motif_id>.pdf         ← if --footprinting (one per motif)
    └── multiome/
        ├── chromvar_scores.csv               ← if --chromvar
        └── diff_motif_activity_top_motifs_<comparison>.csv   ← if --chromvar + --differential_motif_activity
```

---

## File descriptions

### Always present

| File | Written by | Description |
|------|-----------|-------------|
| `calc_scores.log` | all pipeline steps | Full log of the pipeline run with timestamps and runtimes. |

---

### Python scores pipeline (`calc_python_scores/calc_scores.py`)

| File | Condition | Description |
|------|-----------|-------------|
| `adata_st_scores.h5ad` | always | AnnData object for the spatial transcriptomics data. Contains all computed scores (network scores, Squidpy metrics, etc.) in `obsm` and `uns`. Updated later by `add_to_adata.py` and `calc_multiome_scores.py`. |
| `adata_tg_scores.h5ad` | if Tangram was used | AnnData for the Tangram-projected data (single-cell expression projected onto spatial spots). Carries the same score fields as `adata_st_scores.h5ad` plus Tangram-specific data. |
| `adata_map.h5ad` | if Tangram was used | The raw Tangram mapping AnnData object. Rows = dissociated cells, columns = spatial spots; `X` contains the cell-to-spot probability matrix. |
| `adata_map.X.csv` | if Tangram was used | Tangram mapping matrix as CSV (cells × spots). Row names are cell barcodes; column names are spatial spot barcodes. Used as input by the R multiome pipeline for chromVAR projection and footprinting. |
| `adata_map.var.csv` | if Tangram was used | Metadata for spatial spots extracted from the Tangram mapping AnnData (`adata_map.var`). Contains columns such as cluster labels (`leiden`, `cell_type`). Used as spot metadata by the R pipeline. |

---

### R multiome pipeline (`calc_multiome_scores/calc_multiome_scores_test.R`)

Invoked for any job that has `multiome: { use: true }` in the request payload.

| File | Condition | Description |
|------|-----------|-------------|
| `global_motif_analysis.rds` | always (multiome) | The complete Seurat multiome object (original data, not mapped) after all analysis steps. Contains the chromVAR assay, motif information in `peaks@motifs`, and intermediate results in `@misc` (differentially accessible peaks, enriched motifs, differential motif activity). |
| `motif_to_tf.csv` | if `--chromvar`, `--differential_motif_activity`, `--motif_enrichment`, or `--footprinting` | Two-column CSV mapping JASPAR 2024 motif IDs (`motif_id`) to transcription factor names (`TF`). Generated from the JASPAR 2024 database, human CORE collection. |
| `spot_obj_chromvar.rds` | if `--chromvar` + `--differential_motif_activity` | Spot-level Seurat object. The chromVAR deviation scores from dissociated cells are projected onto spatial spots via the Tangram mapping matrix (weighted mean). Contains: `chromvar` assay (motifs × spots), spot metadata from `adata_map.var.csv`, and differential motif activity results in `@misc$diff_motif_activity` (one table per cluster comparison). |
| `multiome/chromvar_scores.csv` | if `--chromvar` | chromVAR deviation scores table. Rows = dissociated cells, columns = JASPAR motif IDs. Produced by `Signac::RunChromVAR`. Read by `add_to_adata.py` to project scores onto spatial spots in the AnnData object. |
| `multiome/diff_motif_activity_top_motifs_<comparison>.csv` | if `--chromvar` + `--differential_motif_activity` | One CSV per cluster comparison (e.g. `0_vs_rest`, `1_vs_rest`). Columns: `motif`, `p_val`, `avg_diff`, `pct.1`, `pct.2`, `p_val_adj`. Contains the top differentially active motifs for that cluster comparison, sorted by effect size. Read by `add_to_adata.py` and stored in `adata.uns["diff_motif_activity_top_motifs"]`. |
| `spot_obj_footprints.rds` | if `--footprinting` | Spot-level Seurat object with spatial footprint profiles stored in `peaks@positionEnrichment`. Each entry has key `footprint_<motif_id>` and contains a (n_spots + 2) × n_positions sparse matrix (positions relative to motif center; last two rows are `expected` and `motif` bias rows for `Signac::PlotFootprint`). |
| `dissociated_obj_footprints.rds` | if `--footprinting` | The dissociated Seurat object after Tn5 insertion bias computation (stored in `peaks@bias`) and Signac `Footprint()`. Can be reused to plot footprints without recomputing them. |
| `footprint_<motif_id>.pdf` | if `--footprinting` (one per motif) | Footprint plot produced by `Signac::PlotFootprint`, showing the spatially-projected insertion enrichment profile of a motif per spatial cluster. Footprints are projections of dissociated-cell ATAC data and are **not** computed from true spatial ATAC fragments. |

---

### `add_to_adata.py` (with `-multiome`)

This script does **not** create new files; it reads the `multiome/` directory and updates the h5ad AnnData files in-place.

| AnnData key added | Source file | Description |
|-------------------|-------------|-------------|
| `obsm["chromvar_spot_scores"]` | `multiome/chromvar_scores.csv` | Spot-level chromVAR scores projected via the Tangram mapping matrix M (spatial spots × motifs; `M.T @ S`). |
| `uns["chromvar_motifs"]` | `multiome/chromvar_scores.csv` | Ordered list of JASPAR motif IDs corresponding to columns in `chromvar_spot_scores`. |
| `uns["diff_motif_activity_top_motifs"][<comparison>]` | `multiome/diff_motif_activity_top_motifs_<comparison>.csv` | Per-comparison DataFrames of top differentially active motifs. |

---

### `calc_multiome_scores.py` (spatial autocorrelation)

This script reads `adata_tg_scores.h5ad`, computes spatial autocorrelation on the projected chromVAR scores, and writes changes back to both `adata_tg_scores.h5ad` and `adata_st_scores.h5ad`.

| AnnData key added | Condition | Description |
|-------------------|-----------|-------------|
| `uns["chromvar_moranI"]` | if `-moranI` | DataFrame of Moran's I spatial autocorrelation statistics for each motif. Columns: `TF`, `motif_id`, `I`, `pval_norm`, `var_norm`, `pval_z_sim`, `pval_sim`, `var_sim`, and FDR-corrected p-values. |
| `uns["chromvar_gearyC"]` | if `-gearyC` | DataFrame of Geary's C spatial autocorrelation statistics for each motif. Columns: `TF`, `motif_id`, `C`, `pval_norm`, `var_norm`, `pval_z_sim`, `pval_sim`, `var_sim`, and FDR-corrected p-values. |

---

## Pipeline execution order

```
1. calc_python_scores/calc_scores.py
      → adata_st_scores.h5ad, adata_tg_scores.h5ad
      → adata_map.h5ad, adata_map.X.csv, adata_map.var.csv   (if Tangram)

2. calc_multiome_scores/calc_multiome_scores_test.R           (if multiome flags set)
      → global_motif_analysis.rds
      → motif_to_tf.csv
      → multiome/chromvar_scores.csv
      → multiome/diff_motif_activity_top_motifs_*.csv
      → spot_obj_chromvar.rds
      → spot_obj_footprints.rds, dissociated_obj_footprints.rds, footprint_*.pdf

3. calc_python_scores/add_to_adata.py  -multiome              (if multiome flags set)
      → updates adata_tg_scores.h5ad / adata_st_scores.h5ad in-place

4. calc_multiome_scores/calc_multiome_scores.py               (if moranI / gearyC flags set)
      → updates adata_tg_scores.h5ad / adata_st_scores.h5ad in-place

5. calc_R_scores/calc_scores.R                                (if R score flags set)
      → Rscores_tg/, Rscores_st/ (temporary, deleted afterwards)

6. calc_python_scores/add_to_adata.py  -Rscores               (if R score flags set)
      → updates adata_tg_scores.h5ad / adata_st_scores.h5ad in-place
```
