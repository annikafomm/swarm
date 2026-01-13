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
    Setzt Zufallssamen für NumPy, random und PYTHONHASHSEED.
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
    Führt PCA, Nachbarn und Louvain-Clustering mit stlearn aus.
    """
    adata.raw = adata

    st.em.run_pca(adata, n_comps=n_comps, random_state=random_state)
    st.pp.neighbors(
        adata,
        n_neighbors=n_neighbors,
        use_rep="X_pca",
        random_state=random_state,
    )
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
    Wählt eine Grid-Kantenlänge n (n x n),
    so dass ~target_cells_per_spot Zellen pro Spot liegen.
    """
    n_cells = adata.n_obs
    n_spots = n_cells / float(target_cells_per_spot)
    n_side = int(np.sqrt(n_spots))

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
    Macht aus Xenium-Zelldaten ein Spot-Level-Grid mit stlearn.
    """
    _set_seed(seed)

    # Wenn kein Cluster-Label vorhanden: Louvain rechnen
    if cluster_key not in adata.obs:
        _run_louvain_clustering(
            adata,
            n_comps=50,
            n_neighbors=25,
            random_state=seed,
            key_added=cluster_key,
        )


    # Grid erzeugen
    grid = st.tl.cci.grid(
        adata,
        n_row=n_spots_side,
        n_col=n_spots_side,
        use_label=cluster_key,
    )

    return grid
