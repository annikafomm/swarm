import scanpy as sc
import numpy as np
from scipy.stats import median_abs_deviation


# filtering functions 

def st_small_filtering(adata):
    # mitochondrial genes
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    adata.obs["mt_outlier"] = adata.obs["pct_counts_mt"] > 20 
    adata = adata[~adata.obs.mt_outlier].copy()

    sc.pp.filter_cells(adata, min_counts=10) # GitHub CoPilot: For spatial transcriptomics (like Xenium) typical values are 10–100.
    sc.pp.filter_genes(adata, min_cells=1)
    

def normalize(adata):
    sc.pp.normalize_total(adata, target_sum=1e4, inplace=True) # Normalize counts per cell
    sc.pp.log1p(adata) # Logarithmize
    sc.pp.pca(adata) # do principal component analysis


def clustering(adata):
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, flavor="igraph", n_iterations=2, directed=False)



def clean_qc_columns(adata):
    """Remove common QC-related columns from AnnData .obs, .var, and .uns."""
    # Columns to remove from obs
    obs_cols_to_remove = [
        'n_genes_by_counts', 'log1p_n_genes_by_counts', 'total_counts', 
        'log1p_total_counts', 'pct_counts_in_top_20_genes', 'total_counts_mt', 
        'log1p_total_counts_mt', 'pct_counts_mt', 'outlier', 'mt_outlier'
    ]
    adata.obs.drop(columns=obs_cols_to_remove, inplace=True, errors="ignore")

    # Columns to remove from var
    var_cols_to_remove = [
        'mt', 'n_cells_by_counts', 'mean_counts', 'log1p_mean_counts', 
        'pct_dropout_by_counts', 'total_counts', 'log1p_total_counts'
    ]
    adata.var.drop(columns=var_cols_to_remove, inplace=True, errors="ignore")

    # Keys to remove from uns
    uns_keys_to_remove = ['outlier_colors', 'mt_outlier_colors', 'log1p']
    for key in uns_keys_to_remove:
        adata.uns.pop(key, None)
