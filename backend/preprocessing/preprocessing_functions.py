import scanpy as sc

# helper functions
def is_outlier(adata, metric: str, nmads: int):
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (np.median(M) + nmads * median_abs_deviation(M) < M)
    return outlier




# filtering functions 

def st_small_filtering(adata):
    # mitochondrial genes
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    adata.obs["mt_outlier"] = adata.obs["pct_counts_mt"] > 20 
    adata = adata[~adata.obs.mt_outlier].copy()

    sc.pp.filter_cells(adata, min_counts=10) # GitHub CoPilot: For spatial transcriptomics (like Xenium) typical values are 10–100.
    sc.pp.filter_genes(adata, min_cells=2) # siehe sc best practices
    #sc.pp.filter_cells(adata, min_genes=...)
    
    


def sc_filtering(adata):
    adata.var_names_make_unique()

    pass



# normalization function

def normalize(adata):
    sc.pp.normalize_total(adata, target_sum=1e4, inplace=True) # Normalize counts per cell
    sc.pp.log1p(adata) # Logarithmize
    sc.pp.pca(adata) # do principal component analysis