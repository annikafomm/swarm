# Pipeline overview:
# 1) Load spatial AnnData (Visium or Xenium).
# 2) Optional preprocessing (filter / normalize).
# 3) If Xenium: build a grid-level AnnData for scoring + keep cell-level AnnData for visualization.
# 4) Optional Tangram -> compute scores on Tangram output.
# 5) Compute scores on spatial data (Visium: original; Xenium: grid-level).
# 6) If Xenium: broadcast grid-level scores back to cells and write cell-level output.


import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning, module="xarray_schema")

import argparse
import os
import scanpy as sc
import squidpy as sq
from scipy import io, sparse
import pandas as pd
import time
import subprocess

import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


from preprocessing.preprocessing_functions import *
from calc_liana import run_liana
from xenium.gridding_mapping import (
    map_cells_to_grid,
    broadcast_grid_to_cells,
)


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

def compute_spatial_scores(adata, description, args, logfile):
    # Calculate spatial scores
    if args.liana or args.centrality_scores or args.co_occurrence or args.nhood_enrichment or args.moranI or args.gearyC:

        if description == "tg":
            log_message("Preparing score calculation for the Tangram output ...", logfile)
        elif description == "st":
            log_message("Preparing score calculation for the Spatial data ...", logfile)

        t0 = time.time()
        sq.gr.spatial_neighbors(adata, coord_type="generic", delaunay=True)
        log_message(f"Spatial neighbors calculated in {format_runtime(t0)}", logfile, 2)

        # liana
        if args.liana:
            if description == "tg":
                cell_prop_key = "tangram_ct_pred"
            elif description == "st":
                cell_prop_key = args.cell_comp_key

            if not (args.grn is None or os.path.exists(args.grn)):
                log_message(f"The GRN file {args.grn} does not exist.", logfile)
            elif not (args.pathway_net is None or os.path.exists(args.pathway_net)):
                log_message(f"The pathway net file {args.pathway_net} does not exist.", logfile)
            else:
                if cell_prop_key not in adata.obsm.keys():
                    log_message(f"'{cell_prop_key}' is not a column in adata.obsm. Please provide a valid cell composition key.", logfile)

                t0 = time.time()
                run_liana(adata, args.grn, args.pathway_net, cell_prop_key)
                log_message(f"LIANA+ scores calculated in {format_runtime(t0)}", logfile, 2)

        # squidpy

        # TODO: Grenze für n_perms

        # check if the cluster key exists in adata.obs if needed
        if args.centrality_scores or args.co_occurrence or args.nhood_enrichment:
            if "leiden" not in adata.obs.keys() and (args.cluster_cs == "leiden" or args.cluster_co == "leiden" or args.cluster_nhood == "leiden"):
                t0 = time.time()
                # neighbors, umap, leiden
                clustering(adata)  # not user configurable, because makeshift solution for when no cluster key is provided
                log_message(f"Clusters calculated in {format_runtime(t0)}", logfile, 2)

        # Compute centrality scores
        if args.centrality_scores:
            if args.cluster_cs not in adata.obs.keys():
                log_message(f"'{args.cluster_cs}' is not a column in adata.obs. Please provide a valid cluster column.", logfile)
            else:
                t0 = time.time()
                sq.gr.centrality_scores(adata, cluster_key=args.cluster_cs, show_progress_bar=True)
                log_message(f"Centrality scores for cluster {args.cluster_cs} calculated in {format_runtime(t0)}", logfile, 2)

        # Compute co-occurrence probability
        if args.co_occurrence:
            if args.cluster_co not in adata.obs.keys():
                log_message(f"'{args.cluster_co}' is not a column in adata.obs. Please provide a valid cluster column.", logfile)
            else:
                t0 = time.time()
                sq.gr.co_occurrence(adata, cluster_key=args.cluster_co, interval = args.interval, n_splits = args.n_splits, show_progress_bar=True)
                log_message(f"Co-occurrence probabilities for cluster {args.cluster_co} calculated in {format_runtime(t0)}", logfile, 2)

        # Compute neighborhood enrichment
        if args.nhood_enrichment:
            if args.cluster_nhood not in adata.obs.keys():
                log_message(f"'{args.cluster_nhood}' is not a column in adata.obs. Please provide a valid cluster column.", logfile)
            elif args.library_key != None and args.library_key not in adata.obs.keys():
                log_message(f"'{args.library_key}' is not a column in adata.obs. Please provide a valid library key.", logfile)
            else:
                t0 = time.time()
                sq.gr.nhood_enrichment(adata, cluster_key=args.cluster_nhood, library_key = args.library_key, seed=42, n_perms=args.n_perms_nhood, show_progress_bar=True)
                log_message(f"Neighborhood enrichment for cluster {args.cluster_nhood} calculated in {format_runtime(t0)}", logfile, 2)

        # Compute Moran's I
        if args.moranI:
            t0 = time.time()
            sq.gr.spatial_autocorr(adata, mode="moran", seed=42, n_perms=args.n_perms_autocorr_mI, transformation=args.n_perms_autocorr_mI is None, two_tailed = args.two_tailed_mI, corr_method = args.corr_method_mI, show_progress_bar=True)
            log_message(f"Moran's I scores calculated in {format_runtime(t0)}", logfile, 2)

        # Compute Geary's C
        if args.gearyC:
            t0 = time.time()
            sq.gr.spatial_autocorr(adata, mode="geary", seed=42, n_perms=args.n_perms_autocorr_gC, transformation=args.n_perms_autocorr_gC is None, two_tailed = args.two_tailed_gC, corr_method = args.corr_method_gC, show_progress_bar=True)
            log_message(f"Geary's C scores calculated in {format_runtime(t0)}", logfile, 2)

        # save AnnData object in file
        log_message("Saving calculations ...", logfile, 2)

    #  TODO: tidy up andata --> delete entries, that are not used further

    t0 = time.time()
    filename = os.path.basename(args.input).replace(".h5ad", f"_{description}_scores.h5ad")
    adata.write(os.path.join(args.outdir, filename))
    log_message(f"AnnData object written in {format_runtime(t0)}", logfile, 4)

    # R scores should be calculated
    if args.R_scores:
        folder_path = os.path.join(args.outdir, f"expr_info_{description}")
        t0 = time.time()

        if not os.path.exists(folder_path):
            os.makedirs(folder_path)
        # Write matrix
        io.mmwrite(os.path.join(folder_path, "expr.mtx"), sparse.csr_matrix(adata.X))
        # Save row names (cells)
        pd.Series(adata.obs_names).to_csv(os.path.join(folder_path, "cells.txt"), index=False, header=False)
        # Save var object
        adata.var.to_csv(os.path.join(folder_path, "var.csv"))

        log_message(f"Expression matrix written in {format_runtime(t0)}", logfile, 4)


def main():
    # Parser
    parser = argparse.ArgumentParser(description="Calculate spatial scores for spatial omics data.")

    # input and output file paths
    parser.add_argument('-input', type=str, required=True, help='Input AnnData file path')
    parser.add_argument('-outdir', type=str, required=True, help='Output dir file path')
    parser.add_argument('-log', type=str, required=True, help='Path to the log file')
    parser.add_argument('-dataset', type=str, choices=['visium', 'xenium'], required=True, help='Type of spatial data')

    # preprocessing options
    parser.add_argument('-filter_st', action='store_true', help='Apply filtering for ST data')
    parser.add_argument('-filter_sc', action='store_true', help='Apply filtering for SC data')
    parser.add_argument('-normalize_st', action='store_true', help='Apply normalization to the ST data')
    parser.add_argument('-normalize_sc', action='store_true', help='Apply normalization to the SC data')

    # tangram
    parser.add_argument("-tangram", action='store_true', help='Compute Tangram')
    parser.add_argument('-sc_path', type=str, help="Path to the single-cell .h5ad file.")
    parser.add_argument('-gene_selection', type=str, choices=['ctg', 'hvg', 'spapros', 'svg', 'None'], default=None, help="Gene selection strategy. Default: use all overlapping genes.")
    parser.add_argument('-cell_label', type=str, default='cell_type', help="Column in adata_sc.obs with cluster/cell annotations (e.g. 'cell_type' or 'cell_subclass').")
    parser.add_argument('-ensembl_col', type=str, default='', help="Column in adata.var with ensembl ids")
    parser.add_argument('-feature_col', type=str, default='', help="Column in adata.var with type of gene")

    # liana
    parser.add_argument("-liana", action='store_true', help='Compute Liana')
    parser.add_argument("-grn", type=str, default=None)
    parser.add_argument("-cell_comp_key", type=str, default="tangram_ct_pred")
    parser.add_argument("-pathway_net", type=str, default=None)

    # squidpy scores
    parser.add_argument('-moranI', action='store_true', help='Compute Moran\'s I')
    parser.add_argument('-gearyC', action='store_true', help='Compute Geary\'s C')
    parser.add_argument('-centrality_scores', action='store_true', help='Compute all centrality scores (closeness, average clustering, degree)')
    parser.add_argument('-co_occurrence', action='store_true', help='Compute co-occurrence probability')
    parser.add_argument('-nhood_enrichment', action='store_true', help='Compute neighborhood enrichment')

    # additional options
    ## centrality scores
    parser.add_argument('-cluster_cs', type=str, default='leiden', help='Cluster key in adata.obs (default: leiden)')
    ## co-occurrence
    parser.add_argument('-interval', type=int, default=50, help='Distance interval for co-occurrence computation (default: 50)')
    parser.add_argument('-n_splits', type=int, default=None, help='Number of splits for co-occurrence computation (default: None)')
    parser.add_argument('-cluster_co', type=str, default='leiden', help='Cluster key in adata.obs (default: leiden)')
    ## nhood enrichment
    parser.add_argument('-library_key', type=str, default=None, help='Key for library information in adata.obs (default: None)')
    parser.add_argument('-n_perms_nhood', type=int, default=1000, help='Number of permutations for neighborhood enrichment')
    parser.add_argument('-cluster_nhood', type=str, default='leiden', help='Cluster key in adata.obs (default: leiden)')
    ## Moran's I
    parser.add_argument('-n_perms_autocorr_mI', type=int, default=None, help='Number of permutations for spatial autocorrelation scores')
    parser.add_argument('-two_tailed_mI', action='store_true', help='Use two-tailed test for spatial autocorrelation scores (default: False)')
    parser.add_argument('-corr_method_mI', type=str, choices=['bonferroni', 'sidak', 'holm-sidak', 'holm', 'simes-hochberg', 'hommel', 'fdr_bh', 'fdr_by', 'fdr_tsbh', 'fdr_tsbky'], default='fdr_bh', help='Correction method for spatial autocorrelation scores (default: benjamini-hochberg)')
    # Geary's C
    parser.add_argument('-n_perms_autocorr_gC', type=int, default=None, help='Number of permutations for spatial autocorrelation scores')
    parser.add_argument('-two_tailed_gC', action='store_true', help='Use two-tailed test for spatial autocorrelation scores (default: False)')
    parser.add_argument('-corr_method_gC', type=str, choices=['bonferroni', 'sidak', 'holm-sidak', 'holm', 'simes-hochberg', 'hommel', 'fdr_bh', 'fdr_by', 'fdr_tsbh', 'fdr_tsbky'], default='fdr_bh', help='Correction method for spatial autocorrelation scores (default: benjamini-hochberg)')


    parser.add_argument('-R_scores', action='store_true', help="Shows if the expression matrix needs to be safed for the calculation of additional scores in R")


    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)


    if args.gene_selection == "None":
        args.gene_selection = None


    # Prepare log file
    logfile = args.log
    log_message(f"Python score pipeline started at {time.strftime('%Y-%m-%d %H:%M:%S')}", logfile)


    # Load the data

    # Spatial Transcriptomics
    if not os.path.exists(args.input):
        log_message(f"The spatial data file {args.input} does not exist.", logfile)
        return
    elif not args.input.endswith('.h5ad'):
        log_message(f"The spatial data file {args.input} has an unsupported file format. Please provide a .h5ad file.", logfile)
        return
    else:
        log_message("Loading ST AnnData object ...", logfile)

        t0 = time.time()
        adata = sc.read_h5ad(args.input)
        log_message(f"AnnData object loaded in {format_runtime(t0)}", logfile, 2)

        # Data objects:
        # - adata: original spatial AnnData (Visium spots OR Xenium cells)
        # - adata_work: the object we compute scores on (Visium: adata, Xenium: grid-level)
        # - adata_cells: Xenium-only, cell-level object kept for visualization + final output
        # - adata_grid: Xenium-only, grid-level object used for scoring

        adata_cells = None          # will only be used for xenium
        adata_grid = None           # will only be used for xenium
        adata_work = adata          # this is what we compute scores on


        st_preprocessed = args.input
        # Preprocessing
        if args.filter_st:
            log_message("Filtering the ST data ...", logfile)
            log_message(f"Number of cells before filtering: {adata.n_obs}", logfile, 2)
            t0 = time.time()

            small_filtering(adata, 'st', args.dataset)

            log_message(f"Number of cells after filtering: {adata.n_obs}", logfile, 2)
            log_message(f"ST data filtered in {format_runtime(t0)}", logfile, 2)

        if args.normalize_st:
            log_message("Normalizing the ST data ...", logfile)
            t0 = time.time()
            # TODO: add check if counts are already normalized
            # is_integer = np.all(np.mod(dense_layer, 1) == 0)
            normalize(adata, args.dataset)
            log_message(f"ST data normalized in {format_runtime(t0)}", logfile, 2)

        if args.filter_st or args.normalize_st:
            st_preprocessed = os.path.join(args.outdir, "st_for_scores.h5ad")
            adata.write(st_preprocessed)

        # If Xenium: create grid-level representation for scoring
        if args.dataset == "xenium":
            log_message("Performing gridding for Xenium data ...", logfile)
            t0 = time.time()
            from xenium.gridding_pipeline import choose_grid_n, gridding_xenium

            # 1) Keep cell-level AnnData (visualization + final output later)
            adata_cells = adata.copy()
            cell_file = os.path.basename(args.input).replace(
                ".h5ad",
                "_xenium_cells.h5ad"
            )
            cell_path = os.path.join(args.outdir, cell_file)
            adata_cells.write(cell_path)
            log_message(f"Cell-level Xenium AnnData written to {cell_path}", logfile, 2)

            # 2) Create grid-level for scoring
            n_ = choose_grid_n(adata, target_cells_per_spot=20)
            log_message(f"Chosen grid size: {n_} x {n_} spots", logfile, 2)

            adata_grid = gridding_xenium(adata, n_spots_side=n_)
            log_message(f"Gridding performed in {format_runtime(t0)}", logfile, 2)

            adata_work = adata_grid

            # 3) Create mapping: each cell -> nearest grid spot
            adata_cells = map_cells_to_grid(
                adata_cells=adata_cells,
                adata_grid=adata_grid,
            )
            # 4) Proceed with all further analyses at the spot level
            st_grid = os.path.join(args.outdir, "st_grid.h5ad")
            adata_work.write(st_grid)

            st_preprocessed = st_grid


        if args.tangram:
            log_message("Prepping Tangram calculations ...", logfile)

            if not os.path.exists(args.sc_path):
                log_message(f"The single cell data file {args.sc_path} does not exist.", logfile)
            elif not args.sc_path.endswith('.h5ad'):
                log_message(f"The single cell data file {args.sc_path} has an unsupported file format. Please provide a .h5ad file.", logfile)
            else:
                log_message("Loading SC AnnData object ...", logfile, 2)

                t0 = time.time()
                adata_sc = sc.read_h5ad(args.sc_path)
                log_message(f"AnnData object loaded in {format_runtime(t0)}", logfile, 4)

                if args.cell_label not in adata_sc.obs.keys():
                    log_message(f"'{args.cell_label}' is not a column in adata.obs (SingleCell). Please provide a valid cell label key.", logfile, 2)
                else:
                    sc_preprocessed = args.sc_path
                    # Preprocessing
                    if args.filter_sc:
                        log_message("Filtering the SC data ...", logfile, 2)
                        log_message(f"Number of cells before filtering: {adata_sc.n_obs}", logfile, 4)
                        t0 = time.time()

                        small_filtering(adata_sc, 'sc', args.dataset)

                        log_message(f"Number of cells after filtering: {adata_sc.n_obs}", logfile, 4)
                        log_message(f"SC data filtered in {format_runtime(t0)}", logfile, 4)

                    if args.normalize_sc:
                        log_message("Normalizing the SC data ...", logfile, 2)
                        t0 = time.time()
                        # TODO: add check if counts are already normalized
                        # is_integer = np.all(np.mod(dense_layer, 1) == 0)
                        normalize(adata_sc, args.dataset)
                        log_message(f"SC data normalized in {format_runtime(t0)}", logfile, 4)

                    if args.filter_sc or args.normalize_sc:
                        sc_preprocessed = os.path.join(args.outdir, "sc_for_tangram.h5ad")
                        adata_sc.write(sc_preprocessed)

                    tangram_out = os.path.join(args.outdir, "tangram_output.h5ad")

                    cmd = [
                        "micromamba", "run", "-n", "tangram",
                        "python", "../backend/calc_python_scores/calc_tangram.py",
                        "--sc_path", sc_preprocessed,
                        "--sp_path", st_preprocessed,
                        "--out_path", tangram_out,
                        "--device", "cpu",
                    ]

                    if args.gene_selection is not None:
                        cmd += ["--gene_selection_mode", args.gene_selection]
                    if args.cell_label:
                        cmd += ["--cell_label", args.cell_label]
                    if args.ensembl_col:
                        cmd += ["--ensembl_col", args.ensembl_col]
                    if args.feature_col:
                        cmd += ["--feature_col", args.feature_col]

                    log_message(f"Running Tangram in tangram_env ...", logfile, 2)
                    t0 = time.time()
                    subprocess.run(cmd, check=True)
                    log_message(f"Tangram script executed in {format_runtime(t0)}", logfile, 4)

                    if os.path.exists(tangram_out):
                        adata_tangram = sc.read_h5ad(tangram_out)
                        compute_spatial_scores(adata_tangram, "tg", args, logfile)
                    else:
                        log_message(
                            f"Tangram output file {tangram_out} not found – skipping Tangram scores.",
                            logfile,
                            2,
                        )
        compute_spatial_scores(adata_work, "st", args, logfile)

        if args.dataset == "xenium":

            # Xenium only: scores were computed on the grid-level (adata_work),
            # so we broadcast grid-level obs/obsm back to each cell using the stored mapping.


            adata_cells = broadcast_grid_to_cells(
                adata_cells=adata_cells,
                adata_grid=adata_work,
                prefix="grid_",
                copy_obs=True,
                copy_obsm=True,
                copy_uns=True,
                copy_varm=True,
                overwrite=False,
            )

            adata_cells.write(
                os.path.join(args.outdir, "xenium_cells_with_grid_scores.h5ad")
            )


        log_message(f"Python score pipeline finished at {time.strftime('%Y-%m-%d %H:%M:%S')}\n", logfile)


if __name__ == "__main__":
    main()
