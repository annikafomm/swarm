import scanpy as sc
import pandas as pd
import numpy as np
import scipy.sparse as sp
from scipy.io import mmwrite
import os

inp = "/workspaces/mopitas-mapra/backend/dgea/GSM6592049_M2_prepro.h5ad"
outdir = "/workspaces/mopitas-mapra/backend/dgea/export_for_seurat"

os.makedirs(outdir, exist_ok=True)

adata = sc.read_h5ad(inp)

# Choice of matrix: prefer "counts" layer, then raw, then main X
X = None
if "counts" in adata.layers:
    X = adata.layers["counts"]
    print("Using adata.layers['counts']")
elif adata.raw is not None:
    X = adata.raw.X
    print("Using adata.raw.X")
else:
    X = adata.X
    print("Using adata.X")

# ensure sparse matrix in CSR format for efficient writing
if not sp.issparse(X):
    X = sp.csr_matrix(X)
else:
    X = X.tocsr()

# Seurat expects features x cells, so transpose if needed
X = X.T.tocsc()

mmwrite(os.path.join(outdir, "matrix.mtx"), X)

# features / genes
var = adata.raw.var if (adata.raw is not None and "counts" not in adata.layers) else adata.var
genes = var.index.to_series()
genes.to_csv(os.path.join(outdir, "genes.tsv"), sep="\t", header=False, index=False)

# barcodes / cells
pd.Series(adata.obs_names).to_csv(os.path.join(outdir, "barcodes.tsv"), sep="\t", header=False, index=False)

# cell metadata
obs = adata.obs.copy()

# convert categorical columns to string for better compatibility with Seurat
for col in obs.columns:
    if str(obs[col].dtype) == "category":
        obs[col] = obs[col].astype(str)

obs.to_csv(os.path.join(outdir, "metadata.csv"))

print("Wrote:", outdir)
print("matrix shape (features x cells):", X.shape)
