import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning, module="xarray_schema")

import argparse
import os
import scanpy as sc
import pandas as pd
import numpy as np
import json
import time
import anndata as ad

# Allow writing nullable strings in anndata
ad.settings.allow_write_nullable_strings = True

from calc_scores import log_message, format_runtime


def combine_files(filename, description, args, logfile):
    t0 = time.time()
    adata_path = os.path.join(args.indir, filename)
    adata = sc.read_h5ad(adata_path)
    log_message(f"AnnData object loaded in {format_runtime(t0)}", logfile, 2)

    t0 = time.time()
    if description == "tg":
        scores_path = os.path.join(args.indir, "Rscores_tg")
    elif description == "st":
        scores_path = os.path.join(args.indir, "Rscores_st")

    for filename in os.listdir(scores_path):
        if filename.endswith('.csv'):
            file_path = os.path.join(scores_path, filename)
            df = pd.read_csv(file_path, index_col=0)

            df_name = filename.replace(".csv", "")
            if df_name in adata.obsm.keys():
                log_message(f"The element with name {df_name} in obsm is overwritten.", logfile, 2)
            adata.obsm[df_name] = df.T

        elif filename.endswith('.json'):
            file_path = os.path.join(scores_path, filename)
            with open(file_path, "r") as f:
                data_dict = json.load(f)

            df_name = filename.replace(".json", "")
            if df_name in adata.uns.keys():
                log_message(f"The element with name {df_name} in uns is overwritten.", logfile, 2)
            adata.uns[df_name] = data_dict

    log_message(f"R score files loaded and added to the AnnData object in {format_runtime(t0)}", logfile, 2)

    print(adata)

    # save AnnData object in file
    log_message("Saving AnnData object ...", logfile, 2)
    t0 = time.time()
    adata.write(adata_path)
    log_message(f"AnnData object written in {format_runtime(t0)}", logfile, 4)


def combine_files_multiome(filename, args, logfile, adata_map_path):
    scores_path = os.path.join(args.indir, "multiome")
    if not os.path.isdir(scores_path):
        log_message(f"Scores directory not found: {scores_path}", logfile, 1)
        return

    if not os.path.isfile(adata_map_path):
        log_message(f"Mapping file not found: {adata_map_path}. Skipping multiome score mapping.", logfile, 1)
        return

    t0 = time.time()
    adata_path = os.path.join(args.indir, filename)
    adata = sc.read_h5ad(adata_path)
    adata_map = sc.read_h5ad(adata_map_path)
    log_message(f"AnnData object loaded in {format_runtime(t0)}", logfile, 2)
    t0 = time.time()

    chromvar_found = False
    files_in_dir = os.listdir(scores_path)
    log_message(f"Files in multiome directory: {files_in_dir}", logfile, 2)

    for filename in files_in_dir:
        if filename.lower() == "chromvar_scores.csv":
            file_path = os.path.join(scores_path, filename)
            chromvar_scores = pd.read_csv(file_path, index_col=0)
            # adata is the AnnData object containing the mapping to spatial data (probabilities)
            M = adata_map.X # cells in dissociated scRNA assay x spatial cells probabilities
            S = np.array(chromvar_scores)[:,1:] # scores for dissociated data  cells x motifs

            M = np.asarray(M, dtype=np.float32)
            S = np.asarray(S, dtype=np.float32)

            spot_chromvar = M.T @ S # spatial cells x motifs
            motif_names = chromvar_scores.columns[1:]

            adata.obsm["chromvar_spot_scores"] = spot_chromvar
            adata.uns["chromvar_motifs"] = list(motif_names)
            chromvar_found = True
            log_message(f"chromvar_spot_scores added to adata.obsm", logfile, 2)

        # differential motif activity
        if filename.lower() == "differential_motif_activity.json":
            file_path = os.path.join(scores_path, filename)
            with open(file_path, "r") as f:
                data_dict = json.load(f)

            adata.uns["differential_motif_activity"] = data_dict
            log_message(f"differential_motif_activity added to adata.uns", logfile, 2)

    if not chromvar_found:
        log_message(f"Warning: chromvar_scores.csv not found in {scores_path}", logfile, 1)

    log_message(f"R score files loaded and added to the AnnData object in {format_runtime(t0)}", logfile, 2)

    # save AnnData object in file
    log_message("Saving AnnData object ...", logfile, 2)
    t0 = time.time()
    adata.write(adata_path)
    log_message(f"AnnData object written in {format_runtime(t0)}", logfile, 4)


def main():
    # Parser

    parser = argparse.ArgumentParser(description="Add R scores to adata file.")

    # input file paths
    parser.add_argument('-indir', type=str, required=True, help='Input dir file path')
    parser.add_argument('-log', type=str, required=True, help='Path to the log file')
    parser.add_argument('-Rscores', action='store_true', help='Whether to merge R scores or not')
    parser.add_argument('-multiome', action='store_true', help='Whether to merge multiome scores or not')

    args = parser.parse_args()

    logfile = args.log
    log_message(f"Combining scores pipeline started at {time.strftime('%Y-%m-%d %H:%M:%S')}", logfile)

    # Load the data
    log_message("Loading score data ...", logfile)

    adatas = []
    for filename in os.listdir(args.indir):
        if filename.endswith('st_scores.h5ad') or filename.endswith('tg_scores.h5ad'):
            adatas.append(filename)

    if len(adatas) <= 0:
        log_message(f"The network scores could not be added to the AnnData object, because there are no .h5ad files in {args.indir}.", logfile)
    if len(adatas) > 2:
        log_message(f"The network scores could not be added to the AnnData object, because there are too many .h5ad files in {args.indir}.", logfile)
    else:
        if args.multiome:
            adata_map_path = os.path.join(args.indir, "adata_map.h5ad")

        for filename in adatas:
            if filename.endswith('tg_scores.h5ad'):
                log_message("Combining files for the Tangram output ...", logfile)
                if args.Rscores:
                    combine_files(filename, "tg", args, logfile)
                if args.multiome:
                    combine_files_multiome(filename, args, logfile, adata_map_path)
            else:
                log_message("Combining files for the Spatial data ...", logfile)
                if args.Rscores:
                    combine_files(filename, "st", args, logfile)
                if args.multiome:
                    combine_files_multiome(filename, args, logfile, adata_map_path)

    log_message(f"Combining scores pipeline finished at {time.strftime('%Y-%m-%d %H:%M:%S')}\n", logfile)


if __name__ == "__main__":
    main()
