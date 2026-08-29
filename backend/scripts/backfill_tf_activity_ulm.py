#!/usr/bin/env python3
"""
Backfill decoupler's ULM TF-activity scores into an already-processed dataset.

Why this exists
---------------
The frontend's TF Activity tab (regulatory-tables-panel.component.ts's actionColumns, and
info.component.ts's own dataset-capability listing, which already documents
`obsm["tf_activity_score_ulm"]` as "optional") expects a *single*, unsuffixed
`adata.obsm["tf_activity_score_ulm"]` / `adata.obsm["tf_activity_padj_ulm"]` pair -- unlike the
AUCell/VIPER/GSVA/ssGSEA regulatory scores, which come in `_genie3`/`_sponge` network-suffixed
variants, ULM is treated as its own single TF-activity-scoring method. No script anywhere in
this repo computes it (confirmed: no reference to "ulm" or "decoupler" outside this file), so
any dataset that never had it hits a 404 the moment its TF Activity tab is opened -- confirmed
for the builtin Heart (Multiome) dataset via `GET /obsm/tf_activity_score_ulm/<gene>`.

Method
------
decoupler's Univariate Linear Model (`dc.mt.ulm`) needs a `(source, target, weight)` network
long-format DataFrame -- reusing the same GENIE3 TF->target network already used for this
dataset's other genie3-suffixed regulatory scores (the CSV at `--network`, columns
`regulatoryGene,targetGene,weight`, renamed to decoupler's expected column names) means this
uses the exact same regulon definitions as AUCell/VIPER for the *same* dataset, rather than
introducing a second, inconsistent TF->target universe. `dc.mt.ulm` needs library-size
normalized + log1p'd expression (not raw counts) -- verified via `adata.X` value inspection
before running (see --skip-normalization-check to bypass if you've already confirmed this).
It writes `adata.obsm["score_ulm"]` / `adata.obsm["padj_ulm"]` in place, which this script then
renames to the `tf_activity_`-prefixed keys the frontend actually looks for.

After writing the AnnData, regenerate the GeoJSON so any changed geometry/meta reaches the
frontend (not required just for this obsm table, since /obsm/{table}/{column} reads the .h5ad
directly rather than going through the GeoJSON -- but keeping the two in sync is still the
supported pattern for anything downstream that does read the GeoJSON's meta).

Example
-------
Dry run against the builtin Heart (Multiome) dataset's live-served ("measured") object --
_resolve_adata_path prefers this file for /obsm/{table}/{column} fetches, so the new obsm table
must live here (not on the Tangram-projected adata_tg_scores.h5ad the network/other genie3
scores were originally computed on) to actually be reachable:

    python3 backend/scripts/backfill_tf_activity_ulm.py \\
        --adata backend/data/heart/plasmidpoop/adata_st_scores_measured.h5ad \\
        --network backend/data/heart/plasmidpoop/genie_network_filtered_tg.csv \\
        --dry-run

For real, with a backup of the input:

    python3 backend/scripts/backfill_tf_activity_ulm.py \\
        --adata backend/data/heart/plasmidpoop/adata_st_scores_measured.h5ad \\
        --network backend/data/heart/plasmidpoop/genie_network_filtered_tg.csv \\
        --backup
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

ad.settings.allow_write_nullable_strings = True


def _fmt(seconds: float) -> str:
    return f"{seconds:.1f}s" if seconds < 60 else f"{seconds / 60:.1f}min"


def _log(message: str) -> None:
    print(f"[backfill-tf-activity-ulm] {message}", flush=True)


def load_network(path: str) -> pd.DataFrame:
    net = pd.read_csv(path)
    rename = {}
    if "regulatoryGene" in net.columns:
        rename["regulatoryGene"] = "source"
    if "targetGene" in net.columns:
        rename["targetGene"] = "target"
    net = net.rename(columns=rename)
    missing = {"source", "target"} - set(net.columns)
    if missing:
        raise SystemExit(
            f"network file {path} is missing required column(s) {missing} after renaming "
            f"(columns found: {list(net.columns)})"
        )
    if "weight" not in net.columns:
        _log("network has no 'weight' column, treating all edges as unweighted (weight=1)")
        net["weight"] = 1.0
    return net[["source", "target", "weight"]]


def check_normalization(adata: ad.AnnData) -> None:
    """decoupler's ULM assumes library-size-normalized, log1p'd expression -- raw integer
    counts violate its own stated assumption and silently produce meaningless scores rather
    than an error, so this is a real check rather than a formality."""
    X = adata.X
    sample = X[:200].toarray() if hasattr(X, "toarray") else np.asarray(X[:200])
    nonzero = sample[sample > 0]
    if nonzero.size == 0:
        _log("WARNING: sampled expression matrix is all zero, cannot verify normalization")
        return
    looks_like_raw_counts = np.allclose(nonzero, np.round(nonzero))
    if looks_like_raw_counts:
        raise SystemExit(
            "adata.X looks like raw integer counts, not library-size-normalized + log1p'd "
            "expression -- decoupler's ULM assumes the latter. Pass --skip-normalization-check "
            "if this is actually fine (e.g. a different but still appropriate transform), or "
            "point --adata at an object whose .X is already normalized."
        )
    _log(f"adata.X looks normalized (sample non-zero range {nonzero.min():.3f}-{nonzero.max():.3f})")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--adata", required=True, help="Path to the .h5ad to update.")
    parser.add_argument("--network", required=True, help="Path to the TF->target network CSV (regulatoryGene,targetGene,weight or source,target,weight).")
    parser.add_argument("--out", default=None, help="Where to write. Defaults to updating --adata in place.")
    parser.add_argument("--tmin", type=int, default=3, help="Minimum number of targets per source TF (decoupler default-equivalent).")
    parser.add_argument("--force", action="store_true", help="Recompute even where tf_activity_score_ulm already exists.")
    parser.add_argument("--backup", action="store_true", help="Copy the input to <name>.pre-backfill.h5ad before overwriting in place.")
    parser.add_argument("--skip-normalization-check", action="store_true", help="Skip the raw-counts sanity check on adata.X.")
    parser.add_argument("--dry-run", action="store_true", help="Report what would be computed, then exit without writing.")
    args = parser.parse_args()

    if not os.path.exists(args.adata):
        raise SystemExit(f"not found: {args.adata}")
    if not os.path.exists(args.network):
        raise SystemExit(f"not found: {args.network}")

    out_path = args.out or args.adata
    in_place = os.path.abspath(out_path) == os.path.abspath(args.adata)

    _log(f"reading {args.adata}")
    t0 = time.time()
    adata = sc.read_h5ad(args.adata)
    _log(f"loaded {adata.n_obs} obs x {adata.n_vars} vars in {_fmt(time.time() - t0)}")

    if not args.force and "tf_activity_score_ulm" in adata.obsm and "tf_activity_padj_ulm" in adata.obsm:
        _log("tf_activity_score_ulm / tf_activity_padj_ulm already present, nothing to do (use --force to recompute)")
        return 0

    net = load_network(args.network)
    _log(f"loaded network: {len(net)} edges, {net['source'].nunique()} regulators")

    if not args.skip_normalization_check:
        check_normalization(adata)

    if args.dry_run:
        _log(f"dry run: would compute ULM with tmin={args.tmin}, nothing written")
        return 0

    import decoupler as dc

    t0 = time.time()
    dc.mt.ulm(adata, net, tmin=args.tmin, verbose=True)
    _log(f"ULM computed in {_fmt(time.time() - t0)}")

    if "score_ulm" not in adata.obsm or "padj_ulm" not in adata.obsm:
        raise SystemExit(
            f"dc.mt.ulm did not produce the expected obsm keys (got: {list(adata.obsm.keys())}) "
            "-- nothing written."
        )
    adata.obsm["tf_activity_score_ulm"] = adata.obsm.pop("score_ulm")
    adata.obsm["tf_activity_padj_ulm"] = adata.obsm.pop("padj_ulm")
    n_tfs = adata.obsm["tf_activity_score_ulm"].shape[1]
    _log(f"tf_activity_score_ulm / tf_activity_padj_ulm: {n_tfs} TFs scored across {adata.n_obs} spots")

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
    return 0


if __name__ == "__main__":
    sys.exit(main())
