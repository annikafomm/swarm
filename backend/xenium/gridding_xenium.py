import stlearn as st
import matplotlib.pyplot as plt
import numpy as np
import random
import os
import scanpy as sc
import warnings

st.settings.set_figure_params(dpi=120)

seed = 0
np.random.seed(seed)
random.seed(seed)
os.environ["PYTHONHASHSEED"] = str(seed)

warnings.filterwarnings("ignore")

base_dir = "/workspaces/mopitas-mapra/backend/data/output-XETG00125__0042682__3802__20240426__200549"

adata = st.ReadXenium(
    feature_cell_matrix_file=os.path.join(base_dir, "cell_feature_matrix.h5"),
    cell_summary_file=os.path.join(base_dir, "cells.csv.gz"),
    scale=1,
    spot_diameter_fullres=15,
    experiment_xenium_file=os.path.join(base_dir, "experiment.xenium"),
)

print(adata)

adata.write_h5ad("xenium.h5ad")


st.pp.filter_genes(adata, min_counts=10)
st.pp.filter_cells(adata, min_counts=10)

adata.raw = adata

st.em.run_pca(adata, n_comps=50, random_state=0)
st.pp.neighbors(adata, n_neighbors=25, use_rep="X_pca", random_state=0)
st.tl.clustering.louvain(adata, random_state=0)


print(adata.obs["louvain"].value_counts())

fig, ax = plt.subplots(figsize=(8, 6))
st.pl.cluster_plot(adata, use_label="louvain", image_alpha=0, size=4, figsize=(10, 10), show_plot=False)
plt.savefig("xenium_louvain.png", dpi=150, bbox_inches="tight")
plt.close()
print("Plot gespeichert als xenium_louvain.png")

# ------------------------------------------------------------------
# Choose grid size n
# ------------------------------------------------------------------
def choose_grid_n(adata, target_cells_per_spot=20, min_n=40, max_n=180):
    """
    Wählt automatisch eine Grid-Kantenlänge n (n x n Grid),
    so dass im Mittel ~target_cells_per_spot pro Spot landen.
    """
    n_cells = adata.n_obs
    n_spots = n_cells / float(target_cells_per_spot)
    n_side = int(np.sqrt(n_spots))

    n_side = max(min_n, n_side)
    n_side = min(max_n, n_side)
    return n_side

n_ = choose_grid_n(adata, target_cells_per_spot=20)
print(f"Verwendete Grid-Größe: {n_} x {n_} Spots")

# ------------------------------------------------------------------
# 4) Gridding
# ------------------------------------------------------------------
st.pp.normalize_total(adata, target_sum=1e4)

grid = st.tl.cci.grid(adata, n_row=n_, n_col=n_, use_label="louvain")
print("Grid shape:", grid.shape)

print(grid.obs["n_cells"].describe())
