import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning, module="xarray_schema")

import argparse
import os
import scanpy as sc
import pandas as pd
import json
import time

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



def main():
    # Parser

    parser = argparse.ArgumentParser(description="Add R scores to adata file.")

    # input file paths
    parser.add_argument('-indir', type=str, required=True, help='Input dir file path')
    parser.add_argument('-log', type=str, required=True, help='Path to the log file')

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
        for filename in adatas:
            if filename.endswith('tg_scores.h5ad'):
                log_message("Combining files for the Tangram output ...", logfile)
                combine_files(filename, "tg", args, logfile)
            else:
                log_message("Combining files for the Spatial data ...", logfile)
                combine_files(filename, "st", args, logfile)

    log_message(f"Combining scores pipeline finished at {time.strftime('%Y-%m-%d %H:%M:%S')}\n", logfile)


if __name__ == "__main__":
    main()
    