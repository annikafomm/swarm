import scanpy as sc
import numpy as np
from scipy.stats import median_abs_deviation


# filtering functions

def small_filtering(adata, data_type, dataset):
    if data_type == "st":
        if dataset == "visium":
            adata = adata[adata.obs["in_tissue"] == True, :]
            print(f"Spots after removing off-tissue spots: {adata.n_obs}")
        mt_threshold = 20
        counts_threshold = 10 # AI: For spatial transcriptomics (like Xenium) typical values are 10–100.
    elif data_type == "sc":
        mt_threshold = 15
        counts_threshold = 200 # AI: For single cell data typical values are 500–1000.
        # Wu et al. threshold?
    else:
        raise ValueError("Small filtering: You can either chose 'st' or 'sc' as type of data.")

    # mitochondrial genes
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=[20], log1p=True)
    adata.obs["mt_outlier"] = adata.obs["pct_counts_mt"] > mt_threshold
    adata = adata[~adata.obs.mt_outlier].copy()

    sc.pp.filter_cells(adata, min_counts=counts_threshold)


def normalize(adata, dataset):
    # TODO: normalize differently for Xenium
    if dataset == "visium":
        sc.pp.normalize_total(adata, target_sum=1e4, inplace=True) # Normalize counts per cell
        sc.pp.log1p(adata) # Logarithmize
        sc.pp.pca(adata) # do principal component analysis
    else:
        sc.pp.normalize_total(adata, target_sum=1e4, inplace=True)


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
