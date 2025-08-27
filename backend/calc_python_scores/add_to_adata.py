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
    
    adata_path = ""
    found_h5ad = False
    for filename in os.listdir(args.indir):
        if filename.endswith('.h5ad'):
            t0 = time.time()
            adata_path = os.path.join(args.indir, filename)
            adata = sc.read_h5ad(adata_path)
            found_h5ad = True
            log_message(f"AnnData object loaded in {format_runtime(t0)}", logfile, 2)
 
    if not found_h5ad:
        log_message(f"The network scores could not be added to the AnnData object, because the .h5ad file in {args.indir} went missing.")
    else:
        t0 = time.time()
        scores_path = os.path.join(args.indir, "Rscores")
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
        log_message("Saving AnnData object ...", logfile)
        t0 = time.time()
        adata.write(adata_path)
        log_message(f"AnnData object written in {format_runtime(t0)}", logfile, 2)

    log_message(f"Combining scores pipeline finished at {time.strftime('%Y-%m-%d %H:%M:%S')}\n", logfile)


if __name__ == "__main__":
    main()
    