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
    parser.add_argument('-indir', type=str, required=True, help='Input dir file path')

    args = parser.parse_args()

    if not (os.path.exists(args.indir) and os.path.isdir(args.indir)):
        raise FileNotFoundError(f"Folder {args.indir} does not exist.")
    
    # Load the data
    print("Reading adata file ...")
    
    adata_path = ""
    found_h5ad = False
    for filename in os.listdir(args.indir):
        if filename.endswith('.h5ad'):
            adata_path = os.path.join(args.indir, filename)
            adata = sc.read_h5ad(adata_path)
            found_h5ad = True
 
    if not found_h5ad:
        raise FileNotFoundError(f"There is no .h5ad file in {args.indir}.")


    for filename in os.listdir(args.indir):
        if filename.endswith('.csv'):
            file_path = os.path.join(args.indir, filename)
            df = pd.read_csv(file_path, index_col=0)

            df_name = filename.replace(".csv", "", 1)
            if df_name not in adata.obsm.keys():
                adata.obsm[df_name] = df.T
            else:
                print(f"There is already an element with name {df_name} in obsm.")
        
        elif filename.endswith('.json'):
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
    adata.write(adata_path, overwrite=True)


if __name__ == "__main__":
    main()
    