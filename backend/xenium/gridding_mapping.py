"""
Design:
- Cells are visualized at true single-cell coordinates (adata_cells.obsm['spatial'])
- All scores are computed at grid/spot level (adata_grid)
- This function copies grid-level results back to each cell based on nearest-spot mapping
- Spatial coordinates are NOT overwritten
"""

import numpy as np
import scanpy as sc
from scipy.spatial import cKDTree
import pandas as pd
from scipy import sparse



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
    copy_uns: bool = True,
    copy_varm: bool = True,
    overwrite: bool = False,
    skip_obsm_keys=("spatial",),
    skip_uns_keys=(),
    skip_varm_keys=(),
) -> sc.AnnData:
    """
    Broadcast grid-level results back to cells using adata_cells.obs[spot_col] mapping.

    - OBS:   adata_grid.obs columns -> adata_cells.obs (optional prefix)
    - OBSM:  adata_grid.obsm tables -> adata_cells.obsm (optional prefix; spatial wird geskippt)
    - UNS:   adata_grid.uns -> adata_cells.uns (gleiche Keys)
    - VARM:  adata_grid.varm -> adata_cells.varm (gleiche Keys; nur wenn var_names kompatibel)
    - OBSP:  wird NICHT kopiert
    """

    if spot_col not in adata_cells.obs:
        raise KeyError(f"{spot_col} not found in adata_cells.obs")

    # Normalize indices to string for safe mapping
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
            # Special case: grid_leiden -> Leiden
            if col == "leiden":
                out_col = "leiden"
            else:
                out_col = f"{prefix}{col}" if prefix else col

            if (not overwrite) and (out_col in adata_cells.obs.columns):
                continue

            adata_cells.obs[out_col] = adata_grid.obs.iloc[grid_pos][col].values

    # ---- OBSM ----
    if copy_obsm:
        for key, mat in adata_grid.obsm.items():
            if key in skip_obsm_keys:
                continue
            out_key = f"{prefix}{key}" if prefix else key
            if (not overwrite) and (out_key in adata_cells.obsm):
                continue

            if isinstance(mat, pd.DataFrame):
                df = mat.iloc[grid_pos].copy()
                df.index = adata_cells.obs_names
                adata_cells.obsm[out_key] = df
                continue

            if sparse.issparse(mat):
                adata_cells.obsm[out_key] = mat[grid_pos, :].copy()
                continue

            arr = np.asarray(mat)
            if arr.ndim != 2 or arr.shape[0] != adata_grid.n_obs:
                continue
            adata_cells.obsm[out_key] = arr[grid_pos, :].copy()

    # ---- UNS (same keys) ----
    if copy_uns:
        for k, v in adata_grid.uns.items():
            if k in skip_uns_keys:
                continue
            if (not overwrite) and (k in adata_cells.uns):
                continue
            adata_cells.uns[k] = v

    # ---- VARM (same keys) ----
    if copy_varm:
        # varm is var-aligned -> require identical var_names
        if not adata_cells.var_names.equals(adata_grid.var_names):
            raise ValueError(
                "Cannot copy .varm because adata_cells.var_names != adata_grid.var_names. "
                "Align variables first (same order) or disable copy_varm."
            )
        for k, v in adata_grid.varm.items():
            if k in skip_varm_keys:
                continue
            if (not overwrite) and (k in adata_cells.varm):
                continue
            adata_cells.varm[k] = v

    return adata_cells

