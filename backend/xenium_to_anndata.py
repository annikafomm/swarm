import scanpy as sc
import pandas as pd
import anndata as ad
import os

path = "datasets"
dir = "output-XETG00125__0042682__3802__20240426__200549"

# Load cell-level matrix
adata = sc.read_10x_h5(os.path.join(path, dir, "cell_feature_matrix.h5"))  # or construct it from .parquet

# Load metadata
cells = pd.read_parquet(os.path.join(path, dir, "cells.parquet"))
adata.obs = cells.set_index("cell_id")

# Save as compressed AnnData
adata.write(os.path.join(path, f"{dir.split('_')[0]}_compact.h5ad"), compression="gzip")
