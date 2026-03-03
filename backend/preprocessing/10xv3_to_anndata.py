import squidpy as sq
import pandas as pd
import numpy as np
import anndata as ad
import scanpy as sc
import argparse

pd.options.mode.string_storage = "python"
ad.settings.allow_write_nullable_strings = True

parser = argparse.ArgumentParser(description="Convert 10x Visium data to AnnData format")
parser.add_argument("--input_dir", type=str, required=True, help="Path to the 10x Visium data directory")
parser.add_argument("--output_file", type=str, required=True, help="Path to save the output AnnData file")

args = parser.parse_args()

path = args.input_dir

visium = sq.read.visium(path)

visium.var_names_make_unique()
visium.obs_names_make_unique()


print(visium)

# Write to file
visium.write_h5ad(args.output_file)