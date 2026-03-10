from __future__ import annotations

import os
import random
import warnings

from typing import Optional

import numpy as np
import scanpy as sc
import stlearn as st

warnings.filterwarnings("ignore")


def _set_seed(seed: int = 0) -> None:
    """
    Set random seeds for NumPy, Python's random module, and PYTHONHASHSEED
    to ensure reproducible grid construction and clustering.
    """
    np.random.seed(seed)
    random.seed(seed)
    os.environ["PYTHONHASHSEED"] = str(seed)


def _run_louvain_clustering(
    adata: sc.AnnData,
    n_comps: int = 50,
    n_neighbors: int = 25,
    random_state: int = 0,
    key_added: str = "louvain",
) -> None:
    """
    Run PCA, compute a kNN graph, and perform Louvain clustering using stlearn.

    This function adds the following to 'adata':
    - adata.raw: snapshot of the input expression matrix
    - adata.obsm["X_pca"]: PCA embedding
    - adata.obsp["connectivities"] / ["distances"]: neighbor graph
    - adata.obs[key_added]: cluster labels
    """
    # Keep a snapshot of the original data
    adata.raw = adata.copy()

    # Dimensionality reduction
    st.em.run_pca(adata, n_comps=n_comps, random_state=random_state)

    # Build neighborhood graph on PCA representation
    st.pp.neighbors(
        adata,
        n_neighbors=n_neighbors,
        use_rep="X_pca",
        random_state=random_state,
    )

    # Community detection with Louvain
    st.tl.clustering.louvain(
        adata,
        random_state=random_state,
        key_added=key_added,
    )


def choose_grid_n(
    adata: sc.AnnData,
    target_cells_per_spot: int = 20,
    min_n: int = 40,
    max_n: int = 180,
) -> int:
    """
    Choose a grid side length 'n' (resulting in an n x n grid)
    such that each spot contains approximately
    'target_cells_per_spot' cells on average.

    The value is clamped between 'min_n' and 'max_n' to avoid
    extremely coarse or extremely fine grids.
    """
    n_cells = adata.n_obs
    n_spots = n_cells / float(target_cells_per_spot)
    n_side = int(np.sqrt(n_spots))

    # Clamp to allowed range
    n_side = max(min_n, n_side)
    n_side = min(max_n, n_side)
    return n_side


def gridding_xenium(
    adata: sc.AnnData,
    n_spots_side: int,
    cluster_key: str = "louvain",
    seed: int = 0,
    target_sum: float = 1e4,
) -> sc.AnnData:
    """
    Convert Xenium single-cell AnnData into a spot-level grid AnnData using stlearn.

    Design:
    - Input 'adata' contains true single-cell coordinates and expression.
    - A regular spatial grid is constructed over the tissue.
    - Each grid spot aggregates nearby cells and uses their cluster labels
      to form a spot-level representation.
    - The returned AnnData represents the grid and is used for all
      downstream spatial scoring (Squidpy, LIANA, etc.).
    """
    _set_seed(seed)

    # If no cluster labels exist, run Louvain clustering
    if cluster_key not in adata.obs:
        _run_louvain_clustering(
            adata,
            n_comps=50,
            n_neighbors=25,
            random_state=seed,
            key_added=cluster_key,
        )


    # Build the spatial grid
    grid = st.tl.cci.grid(
        adata,
        n_row=n_spots_side,
        n_col=n_spots_side,
        use_label=cluster_key,
    )

    grid.var = adata.var.copy()
    grid.var_names = adata.var_names.copy()

    return grid
