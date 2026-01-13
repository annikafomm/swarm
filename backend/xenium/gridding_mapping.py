import numpy as np
import scanpy as sc
from scipy.spatial import cKDTree
import pandas as pd



def map_cells_to_grid(
    adata_cells: sc.AnnData,
    adata_grid: sc.AnnData,
    cell_xy_key: str = "spatial",
    grid_xy_key: str = "spatial",
    out_col: str = "grid_spot",
    out_dist_col: str = "grid_dist",
) -> sc.AnnData:
    """
    Map each Xenium cell to its nearest stlearn grid spot.
    """

    cell_xy = np.asarray(adata_cells.obsm[cell_xy_key])[:, :2]
    grid_xy = np.asarray(adata_grid.obsm[grid_xy_key])[:, :2]

    tree = cKDTree(grid_xy)
    dists, idx = tree.query(cell_xy, k=1)

    adata_cells.obs[out_col] = adata_grid.obs_names[idx].astype(str)
    adata_cells.obs[out_dist_col] = dists

    return adata_cells



def broadcast_grid_to_cells(
    adata_cells: sc.AnnData,
    adata_grid: sc.AnnData,
    spot_col: str = "grid_spot",
    prefix: str = "grid_",
    copy_obs: bool = True,
    copy_obsm: bool = True,
    overwrite: bool = False,
) -> sc.AnnData:
    """
    Broadcast grid-level results back to cells using adata_cells.obs[spot_col] mapping.

    - Copies adata_grid.obs columns -> adata_cells.obs (optionally prefixed)
    - Copies adata_grid.obsm tables -> adata_cells.obsm (optionally prefixed)
    - Does NOT overwrite existing keys unless overwrite=True.
    """

    if spot_col not in adata_cells.obs:
        raise KeyError(f"{spot_col} not found in adata_cells.obs")

    # Normalize indices to string for safe .loc
    grid_index = adata_grid.obs_names.astype(str)
    cell_spots = adata_cells.obs[spot_col].astype(str)

    # Map each cell to a row index in the grid
    grid_pos = pd.Index(grid_index).get_indexer(cell_spots.values)
    if (grid_pos < 0).any():
        missing = int((grid_pos < 0).sum())
        raise ValueError(
            f"{missing} cells map to unknown grid spots. "
            f"Check {spot_col} and that grid obs_names match."
        )

    # ---- OBS ----
    if copy_obs:
        for col in adata_grid.obs.columns:
            out_col = f"{prefix}{col}" if prefix else col
            if (not overwrite) and (out_col in adata_cells.obs.columns):
                continue
            vals = adata_grid.obs.iloc[grid_pos][col].values
            adata_cells.obs[out_col] = vals

    # ---- OBSM ----
    if copy_obsm:
        for key, mat in adata_grid.obsm.items():
            out_key = f"{prefix}{key}" if prefix else key
            if (not overwrite) and (out_key in adata_cells.obsm):
                continue

            # Case 1: pandas DataFrame stored in obsm
            if isinstance(mat, pd.DataFrame):
                df = mat.iloc[grid_pos].copy()
                df.index = adata_cells.obs_names
                adata_cells.obsm[out_key] = df
                continue

            # Case 2: numpy / sparse -> array-like
            arr = np.asarray(mat)
            if arr.ndim != 2 or arr.shape[0] != adata_grid.n_obs:
                # Skip weird shapes safely
                continue

            adata_cells.obsm[out_key] = arr[grid_pos, :].copy()

    return adata_cells

