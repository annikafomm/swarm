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
import os
import gc


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


def expand_spot_adata_to_cells(
    adata_scores: sc.AnnData,
    adata_map: sc.AnnData,
    spot_col: str = "grid_spot",
    keep_map_obs=("imagecol", "imagerow", "grid_spot", "grid_dist"),
    keep_map_obsm=("spatial",),
    chunk_size: int = 200,
    tmp_dir: str | None = None,
    float_dtype=np.float32,
    copy_uns: bool = True,
) -> sc.AnnData:

    if spot_col not in adata_map.obs:
        raise KeyError(f"'{spot_col}' not found in adata_map.obs")

    if tmp_dir is not None:
        os.makedirs(tmp_dir, exist_ok=True)

    cell_spots = adata_map.obs[spot_col].astype(str)
    spot_index = pd.Index(adata_scores.obs_names.astype(str))

    grid_pos = spot_index.get_indexer(cell_spots.values)

    if (grid_pos < 0).any():
        missing = np.unique(cell_spots.values[grid_pos < 0])[:10]
        raise ValueError(
            f"Some cells reference spots not present in adata_scores: {missing}"
        )

    n_cells = adata_map.n_obs
    n_vars = adata_scores.n_vars

    # ------------------------------------------------------------------
    # Sparse expansion helper (memory safe)
    # ------------------------------------------------------------------
    def _expand_sparse(mat):
        mat = mat.tocsr()

        data_parts = []
        indices_parts = []
        indptr = [0]

        for src_row in grid_pos:
            start = mat.indptr[src_row]
            end = mat.indptr[src_row + 1]

            data_parts.append(mat.data[start:end])
            indices_parts.append(mat.indices[start:end])

            indptr.append(indptr[-1] + (end - start))

        data = np.concatenate(data_parts) if data_parts else np.array([], dtype=mat.data.dtype)
        indices = np.concatenate(indices_parts) if indices_parts else np.array([], dtype=mat.indices.dtype)
        indptr = np.asarray(indptr, dtype=np.int64)

        return sparse.csr_matrix((data, indices, indptr), shape=(n_cells, mat.shape[1]))

    # ------------------------------------------------------------------
    # Dense expansion helper
    # ------------------------------------------------------------------
    def _expand_dense(arr, out_name):

        arr = np.asarray(arr)

        dtype = float_dtype if np.issubdtype(arr.dtype, np.floating) else arr.dtype

        if tmp_dir is not None:
            path = os.path.join(tmp_dir, f"{out_name}.dat")
            out = np.memmap(path, mode="w+", dtype=dtype, shape=(n_cells, arr.shape[1]))
        else:
            out = np.empty((n_cells, arr.shape[1]), dtype=dtype)

        for start in range(0, n_cells, chunk_size):
            end = min(start + chunk_size, n_cells)
            out[start:end] = arr[grid_pos[start:end]]

        if isinstance(out, np.memmap):
            out.flush()

        return out

    # ------------------------------------------------------------------
    # Matrix dispatcher
    # ------------------------------------------------------------------
    def _expand_matrix(mat, out_name):

        if isinstance(mat, pd.DataFrame):
            arr = mat.to_numpy()
            cols = mat.columns.copy()

            out = _expand_dense(arr, out_name)

            return pd.DataFrame(out, index=adata_map.obs_names, columns=cols)

        if sparse.issparse(mat):
            return _expand_sparse(mat)

        return _expand_dense(mat, out_name)

    # ------------------------------------------------------------------
    # X
    # ------------------------------------------------------------------
    print("Expanding X ...")
    X_new = _expand_matrix(adata_scores.X, "X")

    # ------------------------------------------------------------------
    # obs
    # ------------------------------------------------------------------
    print("Expanding obs ...")

    obs_new = adata_scores.obs.iloc[grid_pos].copy()
    obs_new.index = adata_map.obs_names.copy()

    for col in keep_map_obs:
        if col in adata_map.obs:
            obs_new[col] = adata_map.obs[col].values

    # ------------------------------------------------------------------
    # var
    # ------------------------------------------------------------------
    var_new = adata_scores.var.copy()

    # ------------------------------------------------------------------
    # Create AnnData
    # ------------------------------------------------------------------
    adata_cells = sc.AnnData(
        X=X_new,
        obs=obs_new,
        var=var_new,
    )

    adata_cells.obs_names = adata_map.obs_names.copy()

    # ------------------------------------------------------------------
    # obsm
    # ------------------------------------------------------------------
    print("Expanding obsm ...")

    for key, value in adata_scores.obsm.items():
        print(f"  obsm[{key}]")
        adata_cells.obsm[key] = _expand_matrix(value, f"obsm_{key}")
        gc.collect()

    # overwrite with real cell spatial coordinates
    for key in keep_map_obsm:
        if key in adata_map.obsm:
            print(f"  overwriting obsm[{key}] from adata_map")
            val = adata_map.obsm[key]
            adata_cells.obsm[key] = val.copy() if hasattr(val, "copy") else np.array(val)

    # ------------------------------------------------------------------
    # layers
    # ------------------------------------------------------------------
    print("Expanding layers ...")

    for key, value in adata_scores.layers.items():
        print(f"  layer[{key}]")
        adata_cells.layers[key] = _expand_matrix(value, f"layer_{key}")
        gc.collect()

    # ------------------------------------------------------------------
    # varm
    # ------------------------------------------------------------------
    print("Copying varm ...")

    for key, value in adata_scores.varm.items():
        adata_cells.varm[key] = value.copy() if hasattr(value, "copy") else value

    # ------------------------------------------------------------------
    # uns
    # ------------------------------------------------------------------
    if copy_uns:
        print("Copying uns ...")
        for key, value in adata_scores.uns.items():
            adata_cells.uns[key] = value

    return adata_cells
