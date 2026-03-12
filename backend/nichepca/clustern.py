import scanpy as sc
import nichepca as npc
import scipy.sparse as sp
import numpy as np

adata = sc.read_h5ad("/workspaces/mopitas-mapra/backend/data/ad_merge.h5ad")

npc.wf.nichepca(adata, knn=25, sample_key = "Sample_ID")
sc.pp.neighbors(adata, use_rep="X_npca")
sc.tl.leiden(adata, resolution=0.5, key_added="leiden_niche")

#sc.pl.spatial(adata, color=["leiden"], size=1, library_id="GSM6592049_M2", img_key=None)

adata.write_h5ad("adata_niche_tangram.h5ad")
