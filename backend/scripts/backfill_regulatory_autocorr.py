#!/usr/bin/env python3
"""
Backfill per-regulatory-score Moran's I / Geary's C into an already-processed dataset.

Why this exists
---------------
`visium_to_geojson.py` builds `meta["global_regulatory_moranI_genie3"]` /
`_gearyC_genie3` / `_moranI_sponge` / `_gearyC_sponge` by reading, for each regulatory
score matrix in `adata.obsm` (any key ending `_genie3` or `_sponge`, e.g.
`aucell_scores_genie3`), a corresponding `adata.uns[f"{key}_moranI"]` /
`adata.uns[f"{key}_gearyC"]` DataFrame -- see visium_to_geojson.py's
`global_regulatory_moranI_genie3` block. Those per-score uns entries are only written by
`compute_regulatory_autocorr()` in calc_python_scores/add_to_adata.py, which only runs as
part of the full upload pipeline. A dataset processed before that step existed, or without
it enabled, has real (non-empty) `global_regulatory_scores_genie3`/`_sponge` (the raw
AUCell/VIPER/GSVA/ssGSEA values) but permanently empty Moran's I/Geary's C tables in the
Regulatory Scores tab -- this was confirmed for the builtin BRCA dataset, which has no
config.json (unlike heart/visual_cortex) and predates this step.

Re-running the whole upload pipeline to recover just this is wasteful for the same reason
`backfill_squidpy_scores.py` exists: it would redo every other score family from scratch.
This script computes only the missing regulatory-score autocorrelation in place, reusing
`compute_regulatory_autocorr`'s exact algorithm (copied here rather than imported, so this
stays independent of calc_scores.py's own heavier import chain -- LIANA, Xenium gridding --
which this backfill has no need for):

  * spatial graph: `sq.gr.spatial_neighbors(adata, coord_type="generic", delaunay=True)`
    (reused if already present)
  * per score matrix: wrap it in a throwaway AnnData sharing the real spatial graph, run
    `sq.gr.spatial_autocorr(mode="moran")` and `mode="geary"` (analytic approximation,
    n_perms=None, matching the pipeline default), store the result back on the real adata
    as `uns[f"{key}_moranI"]` / `uns[f"{key}_gearyC"]`.

After writing the AnnData, regenerate the GeoJSON so the new keys reach the frontend:

    python3 backend/visium_to_geojson.py --adata <adata> --outpath <dir>/hexagons.geojson \\
        --data_type visium

Examples
--------
Dry run against the builtin BRCA dataset:

    python3 backend/scripts/backfill_regulatory_autocorr.py \\
        --adata backend/data/brca_visium/adata.h5ad --dry-run

For real, with a backup of the input:

    python3 backend/scripts/backfill_regulatory_autocorr.py \\
        --adata backend/data/brca_visium/adata.h5ad --backup
"""

import argparse
import os
import shutil
import sys
import time

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq

ad.settings.allow_write_nullable_strings = True


def _fmt(seconds: float) -> str:
    return f"{seconds:.1f}s" if seconds < 60 else f"{seconds / 60:.1f}min"


def _log(message: str) -> None:
    print(f"[backfill-regulatory-autocorr] {message}", flush=True)


def find_regulatory_keys(adata: ad.AnnData) -> list:
    return [key for key in adata.obsm.keys() if key.endswith("_genie3") or key.endswith("_sponge")]


def ensure_spatial_graph(adata: ad.AnnData, force: bool = False) -> bool:
    if "spatial" not in adata.obsm:
        _log("obsm['spatial'] is missing -- spatial coordinates are required. Nothing to do.")
        return False
    if not force and "spatial_connectivities" in adata.obsp:
        _log("spatial graph already present, reusing it")
        return True
    t0 = time.time()
    sq.gr.spatial_neighbors(adata, coord_type="generic", delaunay=True)
    _log(f"spatial graph built in {_fmt(time.time() - t0)}")
    return True


def compute_regulatory_autocorr(adata: ad.AnnData, keys: list, force: bool = False) -> None:
    """Mirrors calc_python_scores/add_to_adata.py's compute_regulatory_autocorr exactly."""
    spatial_connectivities = adata.obsp.get("spatial_connectivities", None)
    spatial_distances = adata.obsp.get("spatial_distances", None)

    for key in keys:
        moran_uns_key, geary_uns_key = f"{key}_moranI", f"{key}_gearyC"
        if not force and moran_uns_key in adata.uns and geary_uns_key in adata.uns:
            _log(f"{key}: already present, skipping (use --force to recompute)")
            continue

        try:
            score_df = adata.obsm[key]
            if not isinstance(score_df, pd.DataFrame):
                score_df = pd.DataFrame(score_df, index=adata.obs_names)

            if score_df.shape[1] == 0:
                _log(f"{key}: matrix has no columns, skipping")
                continue

            tmp = ad.AnnData(X=score_df.to_numpy(dtype=np.float32))
            tmp.obs_names = adata.obs_names.copy()
            tmp.var_names = pd.Index(score_df.columns.astype(str))
            tmp.obsp["spatial_connectivities"] = spatial_connectivities.copy()
            if spatial_distances is not None:
                tmp.obsp["spatial_distances"] = spatial_distances.copy()

            t0 = time.time()
            sq.gr.spatial_autocorr(tmp, mode="moran", n_perms=None, corr_method="fdr_bh", show_progress_bar=False)
            moran_df = tmp.uns.get("moranI", None)
            if isinstance(moran_df, pd.DataFrame):
                adata.uns[moran_uns_key] = moran_df

            sq.gr.spatial_autocorr(tmp, mode="geary", n_perms=None, corr_method="fdr_bh", show_progress_bar=False)
            geary_df = tmp.uns.get("gearyC", None)
            if isinstance(geary_df, pd.DataFrame):
                adata.uns[geary_uns_key] = geary_df

            _log(f"{key}: autocorrelation computed in {_fmt(time.time() - t0)} ({score_df.shape[1]} columns)")
        except Exception as exc:
            _log(f"{key}: FAILED -- {type(exc).__name__}: {exc}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--adata", required=True, help="Path to the .h5ad to update.")
    parser.add_argument("--out", default=None, help="Where to write. Defaults to updating --adata in place.")
    parser.add_argument("--force", action="store_true", help="Recompute even where results already exist.")
    parser.add_argument("--backup", action="store_true", help="Copy the input to <name>.pre-backfill.h5ad before overwriting in place.")
    parser.add_argument("--dry-run", action="store_true", help="Report what is missing and would be computed, then exit without writing.")
    args = parser.parse_args()

    if not os.path.exists(args.adata):
        raise SystemExit(f"not found: {args.adata}")

    out_path = args.out or args.adata
    in_place = os.path.abspath(out_path) == os.path.abspath(args.adata)

    _log(f"reading {args.adata}")
    t0 = time.time()
    adata = sc.read_h5ad(args.adata)
    _log(f"loaded {adata.n_obs} obs x {adata.n_vars} vars in {_fmt(time.time() - t0)}")

    keys = find_regulatory_keys(adata)
    if not keys:
        _log("no *_genie3 / *_sponge regulatory score matrices found in obsm, nothing to do")
        return 0

    todo = [k for k in keys if args.force or f"{k}_moranI" not in adata.uns or f"{k}_gearyC" not in adata.uns]

    _log("current state:")
    for key in keys:
        state = "present" if (f"{key}_moranI" in adata.uns and f"{key}_gearyC" in adata.uns) else "MISSING"
        planned = "will compute" if key in todo else "skip"
        _log(f"    {key:40} {state:8} -> {planned}")

    if args.dry_run:
        _log("dry run, nothing written")
        return 0
    if not todo:
        _log("nothing to do (use --force to recompute)")
        return 0

    if not ensure_spatial_graph(adata, force=args.force):
        return 1

    compute_regulatory_autocorr(adata, todo, force=args.force)

    if in_place and args.backup:
        backup = args.adata.replace(".h5ad", ".pre-backfill.h5ad")
        if os.path.exists(backup):
            _log(f"backup already exists, leaving it alone: {backup}")
        else:
            _log(f"backing up to {backup}")
            shutil.copy2(args.adata, backup)

    for value in adata.obsm.values():
        if hasattr(value, "columns"):
            value.columns = value.columns.astype(str)

    _log(f"writing {out_path}")
    t0 = time.time()
    adata.write(out_path)
    _log(f"written in {_fmt(time.time() - t0)}")
    _log("now regenerate the GeoJSON so the frontend sees the new keys:")
    _log(f"    python3 backend/visium_to_geojson.py --adata {out_path} --outpath <dataset_dir>/hexagons.geojson --data_type visium")
    return 0


if __name__ == "__main__":
    sys.exit(main())
