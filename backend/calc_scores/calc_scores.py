import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning, module="xarray_schema")

import argparse
import os
import scanpy as sc
import squidpy as sq
from scipy import io
import pandas as pd

from preprocessing.preprocessing_functions import *

def main():
    # Parser

    parser = argparse.ArgumentParser(description="Calculate spatial scores for spatial omics data.")

    # input and output file paths
    parser.add_argument('-input', type=str, required=True, help='Input AnnData file path')
    parser.add_argument('-output_dir', type=str, required=True, help='Output dir file path')

    # preprocessing options
    parser.add_argument('-filter', action='store_true', help='Apply filtering')
    parser.add_argument('-normalize', action='store_true', help='Apply normalization')

    # possible scores
    parser.add_argument('-moranI', action='store_true', help='Compute Moran\'s I')
    parser.add_argument('-gearyC', action='store_true', help='Compute Geary\'s C')

    parser.add_argument('-centrality_scores', action='store_true', help='Compute all centrality scores (closeness, average clustering, degree)')    
    parser.add_argument('-co_occurrence', action='store_true', help='Compute co-occurrence probability')
    parser.add_argument('-nhood_enrichment', action='store_true', help='Compute neighborhood enrichment')

    # additional options
    ## centrality scores, co-occurrence, and nhood enrichment
    parser.add_argument('-cluster', type=str, default='leiden', help='Cluster key in adata.obs (default: leiden)')
    ## co-occurrence
    parser.add_argument('-interval', type=int, default=50, help='Distance interval for co-occurrence computation (default: 50)')
    parser.add_argument('-n_splits', type=int, default=None, help='Number of splits for co-occurrence computation (default: None)')
    ## nhood enrichment
    parser.add_argument('-library_key', type=str, default=None, help='Key for library information in adata.obs (default: None)')
    parser.add_argument('-n_perms_nhood', type=int, default=1000, help='Number of permutations for neighborhood enrichment')
    ## Moran's I / Geary's C
    parser.add_argument('-n_perms_autocorr', type=int, default=None, help='Number of permutations for spatial autocorrelation scores')
    parser.add_argument('-two_tailed', action='store_true', help='Use two-tailed test for spatial autocorrelation scores (default: False)')
    parser.add_argument('-corr_method', type=str, default='fdr_bh', help='Correction method for spatial autocorrelation scores (default: benjamini-hochberg)')

    args = parser.parse_args()


    # Load the data
    print("Reading input file ...")
    if not os.path.exists(args.input):
        raise FileNotFoundError(f"Input file {args.input} does not exist.")

    if args.input.endswith('.h5ad'):
        adata = sc.read_h5ad(args.input)
    else:
        raise ValueError("Unsupported file format. Please provide a .h5ad or .zarr file.")

    # Preprocessing

    if args.filter:
        print("Filtering ...")
        print("Number of cells and genes before filtering: ", (adata.n_obs, adata.n_vars))
        st_small_filtering(adata)
        print("Number of cells and genes after filtering: ", (adata.n_obs, adata.n_vars))

    if args.normalize:
        print("Normalization ...")
        # TODO: add check if counts are already normalized
        # is_integer = np.all(np.mod(dense_layer, 1) == 0)
        normalize(adata)

    # Calculate spatial scores


    print("Computing the spatial neighbors ...")
    sq.gr.spatial_neighbors(adata, coord_type="generic", delaunay=True)

    if args.library_key != None and args.library_key not in adata.obs.keys():
            raise ValueError(f"Library key '{args.library_key}' not found in adata.obs. Please provide a valid library key.")
    
    if args.corr_method not in ['bonferroni', 'sidak', 'holm-sidak', 'holm', 'simes-hochberg', 'hommel', 'fdr_bh', 'fdr_by', 'fdr_tsbh', 'fdr_tsbky']:
        raise ValueError(f"Correlation method '{args.corr_method}' not available.")
    
    # TODO: Grenze für n_perms

    # check if the cluster key exists in adata.obs if needed
    if args.centrality_scores or args.co_occurrence or args.nhood_enrichment:
        if args.cluster not in adata.obs.keys() and args.cluster != "leiden":
            raise ValueError(f"Cluster key '{args.cluster}' not found in adata.obs. Please provide a valid cluster key.")
        
        # Compute Leiden clusters if not already present
        if args.cluster == "leiden" and "leiden" not in adata.obs.keys():
            print("Computing Leiden clusters ...")
            # neighbors, umap, leiden
            clustering(adata)  # not user configurable, because makeshift solution for when no cluster key is provided

    # Compute centrality scores
    if args.centrality_scores:
        print("Computing centrality scores ...")
        sq.gr.centrality_scores(adata, cluster_key=args.cluster, show_progress_bar=True)

    # Compute co-occurrence probability
    if args.co_occurrence:
        print("Computing co-occurrence probability ...")
        sq.gr.co_occurrence(adata, cluster_key=args.cluster, interval = args.interval, n_splits = args.n_splits, show_progress_bar=True)

    # Compute neighborhood enrichment
    if args.nhood_enrichment:
        print("Computing neighborhood enrichment ...")
        sq.gr.nhood_enrichment(adata, cluster_key=args.cluster, library_key = args.library_key, seed=42, n_perms=args.n_perms_nhood, show_progress_bar=True)
    

    # Compute Moran's I 
    if args.moranI:
        print(args.attr, args.genes)
        print("Computing Moran's I ...")
        sq.gr.spatial_autocorr(adata, mode="moran", seed=42, n_perms=args.n_perms_autocorr, transformation=args.n_perms_autocorr is None, two_tailed = args.two_tailed, corr_method = args.corr_method, show_progress_bar=True)

    # Compute Geary's C
    if args.gearyC:
        print("Computing Geary's C ...")
        sq.gr.spatial_autocorr(adata, mode="geary", seed=42, n_perms=args.n_perms_autocorr, transformation=args.n_perms_autocorr is None, two_tailed = args.two_tailed, corr_method = args.corr_method, show_progress_bar=True)


    #  TODO: tidy up andata --> delete entries, that are not used further
    
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir, exist_ok=True)
    
    # save AnnData object in file
    print("Saving AnnData object ...")
    adata.write(os.path.join(args.output_dir, f"{os.path.basename(args.input)}_scores.h5ad"))

    # TODO
    if True: # R scores should be calculated
        new_dir = os.path.join(args.output_dir, f"{os.path.basename(args.input)}")
        os.makedirs(new_dir, exist_ok=True)

        # Write matrix
        io.mmwrite(os.path.join(new_dir, "expr.mtx"), adata.X)
        # Save row names (cells)
        pd.Series(adata.obs_names).to_csv(os.path.join(new_dir, "cells.txt"), index=False, header=False)
        # Save var object
        adata.var.to_csv(os.path.join(new_dir, "var.csv"))


if __name__ == "__main__":
    main()
    