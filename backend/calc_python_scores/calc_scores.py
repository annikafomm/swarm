import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning, module="xarray_schema")

import argparse
import os
import scanpy as sc
import squidpy as sq
from scipy import io
import pandas as pd
import time

import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


from preprocessing.preprocessing_functions import *

def log_message(msg, logfile, indent=0):
    prefix = " " * indent
    line = f"{prefix}{msg}"
    print(line)  # console
    with open(logfile, "a") as f:
        f.write(line + "\n")  # write once, no extra blank lines

def format_runtime(t0):
    elapsed = time.time() - t0  # seconds since start
    mins = int(elapsed // 60)
    secs = int(round(elapsed % 60))
    return f"{mins} min {secs} sec"


def main():
    # Parser
    parser = argparse.ArgumentParser(description="Calculate spatial scores for spatial omics data.")

    # input and output file paths
    parser.add_argument('-input', type=str, required=True, help='Input AnnData file path')
    parser.add_argument('-outdir', type=str, required=True, help='Output dir file path')
    parser.add_argument('-log', type=str, required=True, help='Path to the log file')

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

    # Prepare log file
    logfile = args.log
    log_message(f"Python score pipeline started at {time.strftime('%Y-%m-%d %H:%M:%S')}\n", logfile)
    

    # Load the data
    if not os.path.exists(args.input):
        raise FileNotFoundError(f"Input file {args.input} does not exist.")

    if args.input.endswith('.h5ad'):
        log_message("Loading AnnData object ...", logfile)

        t0 = time.time()
        adata = sc.read_h5ad(args.input)
        log_message(f"AnnData object loaded in {format_runtime(t0)}", logfile, 2)
    else:
        raise ValueError("Unsupported file format. Please provide a .h5ad or .zarr file.")

    
    # Preprocessing

    if args.filter:
        log_message("Filtering the ST data ...", logfile)
        log_message(f"Number of cells before filtering: {adata.n_obs}", logfile, 2)
        t0 = time.time()

        st_small_filtering(adata)

        log_message(f"Number of cells after filtering: {adata.n_obs}", logfile, 2)
        log_message(f"ST data filtered in {format_runtime(t0)}", logfile, 2)

    if args.normalize:
        log_message("Normalizing the ST data ...", logfile)
        t0 = time.time()
        # TODO: add check if counts are already normalized
        # is_integer = np.all(np.mod(dense_layer, 1) == 0)
        normalize(adata)
        log_message(f"ST data normalized in {format_runtime(t0)}", logfile, 2)

    # Calculate spatial scores
    if args.centrality_scores or args.co_occurrence or args.nhood_enrichment or args.moranI or args.gearyC:

        if args.library_key != None and args.library_key not in adata.obs.keys():
                raise ValueError(f"Library key '{args.library_key}' not found in adata.obs. Please provide a valid library key.")

        if args.corr_method not in ['bonferroni', 'sidak', 'holm-sidak', 'holm', 'simes-hochberg', 'hommel', 'fdr_bh', 'fdr_by', 'fdr_tsbh', 'fdr_tsbky']:
            raise ValueError(f"Correlation method '{args.corr_method}' not available.")
        

        log_message("Prepping score calculation ...", logfile)

        t0 = time.time()
        sq.gr.spatial_neighbors(adata, coord_type="generic", delaunay=True)
        log_message(f"Spatial neighbors calculated in {format_runtime(t0)}", logfile, 2)

        # TODO: Grenze für n_perms

        # check if the cluster key exists in adata.obs if needed
        if args.centrality_scores or args.co_occurrence or args.nhood_enrichment:
            if args.cluster not in adata.obs.keys() and args.cluster != "leiden":
                raise ValueError(f"Cluster key '{args.cluster}' not found in adata.obs. Please provide a valid cluster key.")

            # Compute Leiden clusters if not already present
            if args.cluster == "leiden" and "leiden" not in adata.obs.keys():
                t0 = time.time()
                # neighbors, umap, leiden
                clustering(adata)  # not user configurable, because makeshift solution for when no cluster key is provided
                log_message(f"Clusters calculated in {format_runtime(t0)}", logfile, 2)


        log_message("Calculating squidpy scores ...", logfile)
        
        # Compute centrality scores
        if args.centrality_scores:
            t0 = time.time()
            sq.gr.centrality_scores(adata, cluster_key=args.cluster, show_progress_bar=True)
            log_message(f"Centrality scores calculated in {format_runtime(t0)}", logfile, 2)


        # Compute co-occurrence probability
        if args.co_occurrence:
            t0 = time.time()
            sq.gr.co_occurrence(adata, cluster_key=args.cluster, interval = args.interval, n_splits = args.n_splits, show_progress_bar=True)
            log_message(f"Co-occurrence probabilities calculated in {format_runtime(t0)}", logfile, 2)

        # Compute neighborhood enrichment
        if args.nhood_enrichment:
            t0 = time.time()
            sq.gr.nhood_enrichment(adata, cluster_key=args.cluster, library_key = args.library_key, seed=42, n_perms=args.n_perms_nhood, show_progress_bar=True)
            log_message(f"Neighborhood enrichment calculated in {format_runtime(t0)}", logfile, 2)




        # Compute Moran's I
        if args.moranI:
            t0 = time.time()
            sq.gr.spatial_autocorr(adata, mode="moran", seed=42, n_perms=args.n_perms_autocorr, transformation=args.n_perms_autocorr is None, two_tailed = args.two_tailed, corr_method = args.corr_method, show_progress_bar=True)
            log_message(f"Moran's I scores calculated in {format_runtime(t0)}", logfile, 2)

        # Compute Geary's C
        if args.gearyC:
            t0 = time.time()
            sq.gr.spatial_autocorr(adata, mode="geary", seed=42, n_perms=args.n_perms_autocorr, transformation=args.n_perms_autocorr is None, two_tailed = args.two_tailed, corr_method = args.corr_method, show_progress_bar=True)
            log_message(f"Geary's C scores calculated in {format_runtime(t0)}", logfile, 2)


    #  TODO: tidy up andata --> delete entries, that are not used further

    # save AnnData object in file
    log_message("Saving calculations ...", logfile)

    t0 = time.time()
    filename = os.path.basename(args.input).replace(".h5ad", "_scores.h5ad")
    adata.write(os.path.join(args.outdir, filename))
    log_message(f"AnnData object written in {format_runtime(t0)}", logfile, 2)

    # TODO
    if True: # R scores should be calculated
        folder_path = os.path.join(args.outdir, "expr_info")
        t0 = time.time()
        if not os.path.exists(folder_path):
            os.makedirs(folder_path)

        # Write matrix
        io.mmwrite(os.path.join(folder_path, "expr.mtx"), adata.X)
        # Save row names (cells)
        pd.Series(adata.obs_names).to_csv(os.path.join(folder_path, "cells.txt"), index=False, header=False)
        # Save var object
        adata.var.to_csv(os.path.join(folder_path, "var.csv"))

        log_message(f"Expression matrix written in {format_runtime(t0)}", logfile, 2)
    
    log_message(f"Python score pipeline finished at {time.strftime('%Y-%m-%d %H:%M:%S')}\n", logfile)


if __name__ == "__main__":
    main()
