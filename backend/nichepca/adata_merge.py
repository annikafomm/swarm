import anndata as ad
import pandas as pd
import numpy as np
import scanpy as sc

ad_ge = sc.read_h5ad("/workspaces/mopitas-mapra/backend/data/tangram/ad_ge_cell_type.h5ad")
adata = sc.read_h5ad('/workspaces/mopitas-mapra/backend/data/GSM6592049_M2_prepro.h5ad')

assert (adata.obs_names == ad_ge.obs_names).all()

# 2) obs-Spalten von ad nach ad_ge kopieren
for col in adata.obs.columns:
    if col not in ad_ge.obs.columns:
        ad_ge.obs[col] = adata.obs[col]

# 3) var-Infos übertragen

common_vars = ad_ge.var_names.intersection(adata.var_names)
ad_vars_for_ge = adata.var.loc[common_vars]
ad_ge.var = ad_ge.var.join(ad_vars_for_ge, how="left", rsuffix="_from_ad")

# 4) obsm aus ad rüberholen (z.B. UMAP, PCA, spatial), aber nicht überschreiben,
for key, val in adata.obsm.items():
    if key not in ad_ge.obsm:
        ad_ge.obsm[key] = val

# 5) obsp (Nachbarschaften) rüberholen – Dimensionen passen (1448 x 1448)
for key, val in adata.obsp.items():
    if key not in ad_ge.obsp:
        ad_ge.obsp[key] = val

# 6) uns rüberholen (z.B. pca, neighbors, umap, spatial), ohne vorhandene zu überschreiben
for key, val in adata.uns.items():
    if key not in ad_ge.uns:
        ad_ge.uns[key] = val
if "feature_is_filtered" in ad_ge.var.columns:
    # Bool-Spalte + NaNs -> zu sauberem Bool machen
    ad_ge.var["feature_is_filtered"] = (
        ad_ge.var["feature_is_filtered"]
        .fillna(False)     # fehlende auf False
        .astype(bool)      # garantiert bool
    )

# --- 5) optional: alle object-Spalten ohne klaren Typ in string umwandeln ---
# (manchmal kommen beim Join komische Typ-Mixe rein)
for c in ad_ge.var.columns:
    if ad_ge.var[c].dtype == "object":
        # wenn es nicht schon ein reiner string ist, in string casten
        ad_ge.var[c] = ad_ge.var[c].astype(str)




ad_merge = ad_ge.copy()
ad_merge.write('/workspaces/mopitas-mapra/backend/data/ad_merge.h5ad')
