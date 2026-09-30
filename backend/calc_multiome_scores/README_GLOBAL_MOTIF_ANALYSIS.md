# Global Motif Analysis Pipeline

## Overview

The Global Motif Analysis pipeline is a comprehensive multi-stage analysis workflow designed to analyze transcription factor (TF) motif activity in multiome (scRNA-seq + scATAC-seq) spatial data. This pipeline integrates R-based chromatin analysis with Python-based spatial statistical computations to provide deep insights into TF-gene regulatory networks in spatial contexts.

The pipeline processes multiome data by:
1. Computing chromVAR motif deviation scores from ATAC peaks
2. Mapping dissociated cell signatures to spatial spots (Tangram)
3. Projecting motif activities to spatial locations
4. Computing spatial autocorrelation statistics (Moran's I, Geary's C)
5. Performing differential motif activity analysis
6. Generating footprinting visualizations

---

## Pipeline Architecture

### Stage 1: R-based Multiome Analysis (`calc_multiome_scores_test.R`)

**Input Files:**
- `--multiome_rds`: Seurat object containing preprocessed RNA + ATAC data (integrated, normalized)
- `--fragments_tsv_gz`: Indexed ATAC fragment file (gzip compressed)
- `--fragments_tsv_gz_tbi`: Tabix index for fragments file (.tbi)

**Computational Steps:**
1. **Motif Database Integration**: Loads JASPAR2024 motifs and adds them to the ATAC assay
2. **ChromVAR Scoring**: Computes motif deviation scores per cell using chromatin accessibility
3. **Differential Accessibility**: Identifies peaks with differential accessibility between cell types
4. **Motif Enrichment**: Finds enriched motifs in differentially accessible peaks
5. **Cell Type-Spot Mapping**: Applies Tangram algorithm to map dissociated cell signatures to spatial spots
6. **ChromVAR Projection**: Projects chromVAR motif scores onto spatial locations using Tangram mapping matrix
7. **Footprinting**: Computes TF footprint patterns in accessible chromatin
8. **Differential Motif Activity**: Identifies motifs with differential activity between spatial cell types

**Output Files:**
- `dissociated_obj_footprints.rds`: Seurat object with chromVAR scores and footprints (cell level)
- `spot_obj_footprints.rds`: Seurat object with projected motif scores (spot level)
- `adata_map.X.csv`: Tangram mapping matrix (cells × spots probability matrix)
- Log files with analysis progress

**Key Parameters:**
- `--genome`: Reference genome (e.g., "hg38" for human)
- `--chromvar`: Enable chromVAR motif deviation analysis
- `--footprinting`: Enable footprint computation
- `--differential_motif_activity`: Enable differential motif activity analysis

---

### Stage 2: AnnData Integration (`add_to_adata.py`)

**Purpose:** Bridge R-based Seurat objects and Python-based AnnData format

**Input Requirements:**
- Spot-level Seurat RDS files from Stage 1
- AnnData object from Python pipeline

**Processing:**
- Loads projected motif scores from RDS files
- Integrates scores into AnnData obsm (multi-dimensional observations)
- Ensures metadata consistency across formats
- Creates motif-to-TF mapping files

**Output:**
- Updated AnnData object with chromVAR projections stored in `obsm["chromvar_spot_scores"]`
- Motif annotations in `uns["chromvar_motifs"]`

---

### Stage 3: Spatial Correlation Analysis (`calc_multiome_scores.py`)

**Purpose:** Compute spatial statistics for motif activities

**Analytical Methods:**

1. **Moran's I Autocorrelation** (`-moranI`):
   - Tests spatial autocorrelation of motif activities across spots
   - Output: I statistic, p-values (parametric, permutation-based), FDR correction
   - Maps results to transcription factors via motif-to-TF mapping

2. **Geary's C Autocorrelation** (`-gearyC`):
   - Alternative spatial autocorrelation measure (more sensitive to local clusters)
   - Similar output structure to Moran's I
   - Complementary to Moran's I for validation

**Parameters:**
- `-n_perms_autocorr_mI`: Number of permutations for Moran's I (default: 1000)
- `-n_perms_autocorr_gC`: Number of permutations for Geary's C (default: 1000)
- `-corr_method_mI/gC`: Multiple testing correction method (default: "fdr_bh")
- `-two_tailed_mI/gC`: Enable two-tailed statistical tests

**Output:**
- `chromvar_moranI` in `adata.uns`: Moran's I results with TF annotations
- `chromvar_gearyC` in `adata.uns`: Geary's C results with TF annotations
- Final AnnData objects: `adata_tg_scores.h5ad` and `adata_st_scores.h5ad`

---

## Usage

### Running the Complete Pipeline via `app.py`

The pipeline is orchestrated through the `calculate_scores_helper()` function in `app.py`. When `multiome_params` are detected, the following subprocess chain is executed:

```bash
# Parameters are passed as lists that get concatenated
multiome_params = ["--chromvar", "--footprinting", "--genome", "hg38", ...]
multiome_params_py = ["-chromvar", "-moranI", "-gearyC", ...]

# Stage 1: R-based Analysis
Rscript /workspaces/swarm/backend/calc_multiome_scores/calc_multiome_scores_test.R \
    --outdir {out_dir} \
    --log {log_file} \
    {multiome_params}

# Stage 2: AnnData Integration
python3 /workspaces/swarm/backend/calc_python_scores/add_to_adata.py \
    -indir {out_dir} \
    -log {log_file} \
    -multiome

# Stage 3: Spatial Statistics
python3 /workspaces/swarm/backend/calc_multiome_scores/calc_multiome_scores.py \
    --dir {out_dir} \
    --log {log_file} \
    {multiome_params_py}
```

### Standalone Execution

#### Step 1: Run R-based Multiome Analysis

```bash
Rscript calc_multiome_scores_test.R \
    --outdir /path/to/output \
    --log /path/to/output/calc_scores.log \
    --multiome_rds /path/to/multiome_object.rds \
    --fragments_tsv_gz /path/to/fragments.tsv.gz \
    --fragments_tsv_gz_tbi /path/to/fragments.tsv.gz.tbi \
    --genome hg38 \
    --chromvar \
    --footprinting \
    --differential_motif_activity
```

#### Step 2: Integrate into AnnData

```bash
python3 /workspaces/swarm/backend/calc_python_scores/add_to_adata.py \
    -indir /path/to/output \
    -log /path/to/output/calc_scores.log \
    -multiome
```

#### Step 3: Compute Spatial Correlations

```bash
python3 calc_multiome_scores.py \
    --dir /path/to/output \
    --log /path/to/output/calc_scores.log \
    -chromvar \
    -moranI \
    -gearyC \
    -n_perms_autocorr_mI 1000 \
    -n_perms_autocorr_gC 1000
```

---

## Additional Utilities

### Compute Additional Footprints

For post-hoc computation of footprints for specific motifs without rerunning the entire pipeline:

```bash
Rscript compute_additional_footprints.R \
    --outdir /path/to/output \
    --motifs MA1638.2,MA0006.2 \
    --cluster_by cell_type \
    --save_rds FALSE
```

**Parameters:**
- `--motifs`: Comma-separated JASPAR motif IDs
- `--cluster_by`: Metadata column for spatial grouping
- `--cluster_by_mult`: Multiple clustering types (overrides `--cluster_by`)
- `--save_rds`: Save updated RDS files with footprints

---

## Data Formats and Specifications

### Input Data Requirements

1. **Multiome RDS Object**:
   - Seurat object with RNA and ATAC assays
   - Must be preprocessed: log-normalized RNA, LSI reduced ATAC
   - Metadata: cell type annotations, clustering labels
   - Example: `healthy_breast_preprocessed_RNA_peaks_only.rds`

2. **Fragment File**:
   - TSV format (gzip compressed): `chr`, `start`, `end`, `barcode`, `count`
   - Must be indexed with Tabix (.tbi file)
   - One row per unique barcode, position, fragment count

### Output Data Specifications

1. **Spot-level Motif Scores** (stored in AnnData):
   ```
   adata.obsm["chromvar_spot_scores"]: (n_spots, n_motifs)
   adata.uns["chromvar_motifs"]: List of motif names
   adata.uns["chromvar_moranI"]: DataFrame with columns:
       - TF: Transcription factor name
       - motif_id: JASPAR motif identifier
       - I: Moran's I statistic
       - pval_norm, pval_z_sim, pval_sim: Various p-values
       - pval_*_fdr_bh: FDR-corrected p-values
   adata.uns["chromvar_gearyC"]: Similar structure to moranI
   ```

2. **RDS Objects**:
   - `dissociated_obj_footprints.rds`: (n_motifs, n_cells) chromVAR assay
   - `spot_obj_footprints.rds`: (n_motifs, n_spots) chromVAR assay + spatial metadata

### File Naming Convention

Outputs follow a consistent naming pattern:
- RDS files: `{cell_type}_obj_footprints.rds`
- Mapping matrix: `adata_map.X.csv`
- AnnData objects: `adata_tg_scores.h5ad`, `adata_st_scores.h5ad`
- Logs: `calc_scores.log`

---

## Dependencies

### R Packages
- `Seurat >= 5.0`: Single-cell data manipulation
- `Signac >= 1.12`: ATAC-seq analysis
- `JASPAR2024`: Transcription factor motif database
- `TFBSTools`: Motif scanning and annotation
- `BSgenome.Hsapiens.UCSC.hg38`: Reference genome
- `chromVAR`: Motif deviation scoring
- `dplyr`, `tibble`: Data manipulation
- `patchwork`, `ggplot2`: Visualization

### Python Packages
- `scanpy >= 1.9`: Single-cell data analysis
- `squidpy >= 1.3`: Spatial transcriptomics analysis
- `anndata >= 0.9`: Annotated data format
- `pandas`: Data manipulation
- `scipy`: Scientific computing
- `numpy`: Numerical computing

---

## Performance Considerations

### Computational Requirements

- **Time**: Stage 1 (R analysis) typically takes 2-8 hours depending on data size
- **Memory**: ~64 GB RAM recommended for large datasets
- **Disk**: Temporary files can consume 50-100 GB during processing

### Parallelization

- Footprinting uses `parallel::mclapply()` with auto-detected cores
- Set `FOOTPRINT_WORKERS` environment variable to limit workers:
  ```bash
  export FOOTPRINT_WORKERS=4
  Rscript compute_additional_footprints.R ...
  ```

### Optimization Tips

1. Pre-filter genes/peaks to focus on biologically relevant features
2. Use `--genome` parameter appropriate for your data
3. Adjust permutation numbers (`-n_perms_autocorr_*`) based on computational budget
4. Run spatial statistics with/without two-tailed tests based on hypothesis

---

## Troubleshooting

### Common Issues

| Issue | Cause | Solution |
|-------|-------|----------|
| Motifs not added to ATAC assay | JASPAR database not downloaded | Install JASPAR2024 package: `install.packages("JASPAR2024")` |
| Fragment file errors | Invalid tabix index or format | Regenerate index: `tabix -p gff fragments.tsv.gz` |
| Memory errors during projection | Dataset too large for RAM | Split analysis by cell type batches |
| Missing motif-to-TF mapping | File not generated in add_to_adata.py | Check intermediate output directory for mapping files |

### Debug Mode

Enable verbose logging by inspecting log files:
```bash
tail -f /path/to/output/calc_scores.log
```

---

## Citation and References

This pipeline integrates multiple established methods:
- **ChromVAR**: Schep et al., *Nature Methods* (2017)
- **Tangram**: Biancalani et al., *Science* (2021)
- **Signac**: Stuart, Butler, et al., *Nature Biotechnology* (2021)
- **Squidpy**: Palla et al., *Nature Methods* (2022)
- **Moran's I**: Moran, P. A. (1950)
- **Geary's C**: Geary, R. C. (1954)

---

## Authors and Contact

For questions, issues, or contributions, please contact the bioinformatics team or submit an issue to the repository.

---

## Version History

- **v1.0** (2026-04-12): Initial comprehensive documentation
- Previous: Incremental development of multiome analysis components

---

## License

[Add appropriate license information]
