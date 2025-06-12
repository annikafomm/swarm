"""
Source: sc-best-practices book: Heumos, L., Schaar, A.C., Lance, C. et al. 
        Best practices for single-cell analysis across modalities. Nat Rev Genet (2023). 
        https://doi.org/10.1038/s41576-023-00586-w
"""

# Quality control

import numpy as np
import scanpy as sc
import seaborn as sns
from scipy.stats import median_abs_deviation

sc.settings.verbosity = 0
sc.settings.set_figure_params(
    dpi=80,
    facecolor="white",
    frameon=False,
)

import os
import spatialdata as sd

sdata = sd.read_zarr(os.path.join("..", "datasets", "xenium.zarr"))
adata = sdata.tables["table"]
adata.var_names_make_unique()

## Filtering low quality cells

adata.var["mt"] = adata.var_names.str.startswith("MT-" or "mt-")
print(adata.var["mt"].sum(), "mitochondrial genes")

sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=[20], log1p=True)


def is_outlier(adata, metric: str, nmads: int):
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (np.median(M) + nmads * median_abs_deviation(M) < M)
    return outlier

adata.obs["outlier"] = (
    is_outlier(adata, "log1p_total_counts", 5)
    | is_outlier(adata, "log1p_n_genes_by_counts", 5)
    | is_outlier(adata, "pct_counts_in_top_20_genes", 5)
)

adata.obs["mt_outlier"] = is_outlier(adata, "pct_counts_mt", 3) | (
    adata.obs["pct_counts_mt"] > 8
)


print(f"Total number of cells: {adata.n_obs}")
adata = adata[(~adata.obs.outlier) & (~adata.obs.mt_outlier)].copy()
print(f"Number of cells after filtering of low quality cells: {adata.n_obs}")