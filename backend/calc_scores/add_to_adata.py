import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning, module="xarray_schema")

import argparse
import os
import scanpy as sc
import pandas as pd
import json

def main():
    # Parser

    parser = argparse.ArgumentParser(description="Add R scores to adata file.")

    # input file paths
    parser.add_argument('-adata', type=str, required=True, help='Input AnnData file path')
    parser.add_argument('-indir', type=str, required=True, help='Input dir file path')

    args = parser.parse_args()


    # Load the data
    print("Reading adata file ...")
    if not os.path.exists(args.adata):
        raise FileNotFoundError(f"AnnData file {args.adata} does not exist.")

    if args.adata.endswith('.h5ad'):
        adata = sc.read_h5ad(args.adata)
    else:
        raise ValueError("Unsupported file format. Please provide a .h5ad or .zarr file.")

    
    if not (os.path.exists(args.indir) and os.path.isdir(args.indir)):
        raise FileNotFoundError(f"Folder {args.indir} does not exist.")
    

    for filename in os.listdir(args.indir):
        if filename.endswith('.csv'):
            file_path = os.path.join(args.indir, filename)
            df = pd.read_csv(file_path, index_col=0)

            df_name = filename.replace(".csv", "", 1)
            if df_name not in adata.obsm.keys():
                adata.obsm[df_name] = df.T
            else:
                print(f"There is already an element with name {df_name} in obsm.")
    
    
    for filename in os.listdir(args.indir):
        if filename.endswith('.json'):
            file_path = os.path.join(args.indir, filename)
            with open(file_path, "r") as f:
                data_dict = json.load(f)

            df_name = filename.replace(".json", "", 1)
            if df_name not in adata.uns.keys():
                adata.uns[df_name] = data_dict
            else:
                print(f"There is already an element with name {df_name} in uns.")

    # save AnnData object in file
    print("Saving AnnData object ...")
    adata.write(args.adata, overwrite=True)


if __name__ == "__main__":
    main()
    