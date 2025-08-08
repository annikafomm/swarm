import argparse
import os
import numpy as np
import scanpy as sc
from scipy.stats import median_abs_deviation

from preprocessing_functions import *


def is_outlier(adata, metric: str, nmads: int):
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M)
    return outlier

def main():
    # Parser

    parser = argparse.ArgumentParser(description="Preprocess spatial omics data.")

    # input and output file paths
    parser.add_argument('-input', type=str, required=True, help='Input AnnData file path')
    parser.add_argument('-mt_ensembl_ids', type=str, default="mt_ensembl_ids.txt", help='Path to file containing mitochondrial gene Ensembl IDs')  
    parser.add_argument('-output', type=str, default=None, help='Output AnnData file path')

    parser.add_argument('-mt_col', type=str, default=None, help='Column name with gene names as symbols or Ensembl IDs in adata.var. If not provided, the script will look at var_names.')

    args = parser.parse_args()

    # Load the input data
    if not os.path.exists(args.input):
        raise FileNotFoundError(f"Input file {args.input} does not exist.")
    adata = sc.read_h5ad(args.input)
    
    # Read in set of mt ensembl ids
    mt_ensembl_ids = set()  
    if not os.path.exists(args.mt_ensembl_ids):
        raise FileNotFoundError(f"Input file {args.mt_ensembl_ids} does not exist.")
    with open(args.mt_ensembl_ids, "r") as f:
        mt_ensembl_ids = set(line.strip() for line in f)

    # label mitochondrial genes
    if args.mt_col is None:
        adata.var["mt"] = adata.var_names.str.startswith("MT-") | adata.var_names.str.startswith("mt-")
        if adata.var["mt"].sum() == 0:
            print("No mitochondrial genes found by prefix. Checking if var_names are Ensembl IDs...")
            adata.var["mt"] = adata.var_names.isin(mt_ensembl_ids)
        print(f"Number of mitochondrial genes: {adata.var['mt'].sum()}")
    elif args.mt_col in adata.var.columns:
        adata.var["mt"] = adata.var[args.mt_col].str.startswith("MT-") | adata.var[args.mt_col].str.startswith("mt-")
        if adata.var["mt"].sum() == 0:
            print("No mitochondrial genes found by prefix. Checking if column contains Ensembl IDs...")
            adata.var["mt"] = adata.var[args.mt_col].isin(mt_ensembl_ids)
        print(f"Number of mitochondrial genes: {adata.var['mt'].sum()}")
    else:
        raise ValueError(f"Column {args.mt_col} not found in adata.var. Please provide a valid column name.")

    if adata.var["mt"].sum() == 0:
        print("No mitochondrial genes found. Please check the input file and the mt_ensembl_ids file.")

    # calculate the QC covariates
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=[20],log1p=True)

    # relative quality filtering

    adata.obs["outlier"] = (
        is_outlier(adata, "log1p_total_counts", 5)
        | is_outlier(adata, "log1p_n_genes_by_counts", 5)
        | is_outlier(adata, "pct_counts_in_top_20_genes", 5)
    )
    
    # mitochondrial gene fraction
    adata.obs["mt_outlier"] = is_outlier(adata, "pct_counts_mt", 3) | (
    adata.obs["pct_counts_mt"] > 8
    )
    
    # Statistics on how many cells and genes are filtered out
    print(f"Total number of cells: {adata.n_obs}")
    adata = adata[(~adata.obs.outlier) & (~adata.obs.mt_outlier)].copy()
    print(f"Number of cells after filtering of low quality cells: {adata.n_obs}")

    print(f"Total number of genes: {adata.n_vars}")
    sc.pp.filter_genes(adata, min_cells=2)
    print(f'Number of genes after filtering: {adata.n_vars}')

    # Normalize the data
    normalize(adata)

    # Calculate clustering
    sc.pp.neighbors(adata)
    sc.tl.leiden(adata, flavor="igraph", n_iterations=2, directed=False)

    basename = os.path.splitext(os.path.basename(args.input))[0]

    # Write the output file
    if args.output is None:
        args.output = f"{basename}_prepro.h5ad"
    adata.write(args.output)


    sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt", save=f"_{basename}.png", show=False)

if __name__ == "__main__":
    main()