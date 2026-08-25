#!/usr/bin/env python3
"""
Backfill the Squidpy score family into an already-processed dataset.

Why this exists
---------------
The Squidpy scores are only computed during the original upload run, gated on the
`scores.squidpy` flag. A dataset processed with that flag off (for example the builtin
`Heart (Multiome)`, whose config has `"squidpy": false`) is missing:

  * `obs["leiden"]`                        -> no clustering at all
  * `uns["leiden_centrality_scores"]`      -> Cluster Information centrality block
  * `uns["leiden_co_occurrence"]`          -> Co-occurrence tab + `meta.interval`
  * `uns["leiden_nhood_enrichment"]`       -> Cluster Information neighbourhood chart
  * `uns["moranI"]` / `uns["gearyC"]`      -> Gene Expression sorting, and the gene the
                                              map opens on

Re-running the whole upload pipeline to recover those is wasteful (it would redo Tangram,
chromVAR, footprinting and GRN evaluation), and for a builtin dataset there is no upload to
re-run in the first place. This script computes just the missing pieces in place.

It deliberately mirrors `calc_python_scores/calc_scores.py` rather than inventing its own
parameters, so a backfilled dataset is indistinguishable from one processed with
`scores.squidpy` enabled:

  * spatial graph:  `sq.gr.spatial_neighbors(adata, coord_type="generic", delaunay=True)`
  * clustering:     `sc.pp.neighbors` -> `sc.tl.umap` -> `sc.tl.leiden(flavor="igraph",
                    n_iterations=2, directed=False)`  (same as `preprocessing.clustering()`)
  * scores:         `sq.gr.centrality_scores`, `sq.gr.co_occurrence`,
                    `sq.gr.nhood_enrichment`, `sq.gr.spatial_autocorr`

After writing the AnnData, regenerate the GeoJSON so the new keys reach the frontend:

    python3 visium_to_geojson.py --adata <adata> --outpath <dir>/hexagons.geojson \\
        --data_type visium

Note on the cluster key: `visium_to_geojson.py` only ever reads `leiden_centrality_scores`,
`leiden_co_occurrence` and `leiden_nhood_enrichment` by those exact names. Computing under a
different `--cluster-key` therefore produces `uns["<key>_*"]` entries that the frontend will
never display, unless you also pass `--alias-to-leiden`.

Examples
--------
Backfill the builtin heart dataset (dry run first):

    python3 backend/scripts/backfill_squidpy_scores.py \\
        --adata backend/data/heart/plasmidpoop/adata_tg_scores.h5ad --dry-run

    python3 backend/scripts/backfill_squidpy_scores.py \\
        --adata backend/data/heart/plasmidpoop/adata_tg_scores.h5ad --backup

Group the cluster scores by an existing annotation instead of a fresh Leiden run, and file
them under the names the GeoJSON reads:

    python3 backend/scripts/backfill_squidpy_scores.py --adata <adata> \\
        --cluster-key cell_type --alias-to-leiden
"""

import argparse
import os
import shutil
import sys
import time

import scanpy as sc
import squidpy as sq


def _fmt(seconds: float) -> str:
    return f"{seconds:.1f}s" if seconds < 60 else f"{seconds / 60:.1f}min"


def _log(message: str) -> None:
    print(f"[backfill] {message}", flush=True)


def ensure_spatial_graph(adata, force: bool = False) -> None:
    """Build the Delaunay spatial graph the cluster scores are computed on.

    Matches calc_scores.py exactly. Squidpy stores the result in obsp, with its parameters
    in uns["spatial_neighbors"]; both are needed by the scores below.
    """
    if "spatial" not in adata.obsm:
        raise SystemExit(
            "obsm['spatial'] is missing — spatial coordinates are required for every "
            "Squidpy score. Nothing to do."
        )
    if not force and "spatial_neighbors" in adata.uns and "spatial_connectivities" in adata.obsp:
        _log("spatial graph already present, reusing it")
        return
    t0 = time.time()
    sq.gr.spatial_neighbors(adata, coord_type="generic", delaunay=True)
    _log(f"spatial graph built in {_fmt(time.time() - t0)}")


def ensure_leiden(adata, force: bool = False, n_hvg: int = 0, resolution: float = None) -> None:
    """Compute obs['leiden'].

    With no --n-hvg/--resolution this mirrors preprocessing.clustering() exactly:
    expression-space neighbours (not the spatial graph), then UMAP, then Leiden at scanpy's
    default resolution, with the PCA computed on the fly over *all* genes because
    obsm['X_pca'] is absent on a dataset that never clustered.

    That default is faithful but weak, and the pipeline says so itself — the call is commented
    "makeshift solution for when no cluster key is provided". Two knobs make it defensible:

      * `--n-hvg N` restricts the PCA to the N most variable genes. Running it over all ~27k
        genes lets uninformative low-variance genes dominate the distance metric.
      * `--resolution R` controls how finely Leiden cuts the graph. The default of 1.0
        over-partitions when the tissue has only a handful of real states: on the heart
        builtin it produced 20 near-equal clusters (313, 239, 225, 223, ... 66, 54) against
        6 expert niches, with only 70.8% mean purity. A flat size distribution like that is
        the signature of cutting a *continuum* rather than finding communities.

    Note that clustering Tangram-projected expression is inherently harder: projection is a
    mapping-weighted average, so it smooths neighbouring spots together and weakens exactly
    the separation Leiden looks for. Prefer running this on the measured spatial object.
    """
    if not force and "leiden" in adata.obs.columns:
        _log("obs['leiden'] already present, keeping it")
        return

    # A dataset that never clustered can still carry a stale, incomplete uns["neighbors"]
    # (the heart builtin has one holding only {"params": {"method": ...}}, with an empty
    # obsp). scanpy's Neighbors.__init__ treats the presence of that key as "a graph already
    # exists" and dies on the missing params.n_neighbors, so clear it out first rather than
    # letting the run fail on stale metadata.
    stale_uns = "neighbors" in adata.uns and "n_neighbors" not in (
        adata.uns["neighbors"].get("params", {}) or {}
    )
    missing_graph = not {"distances", "connectivities"} & set(adata.obsp.keys())
    if "neighbors" in adata.uns and (stale_uns or missing_graph):
        _log("dropping stale/incomplete uns['neighbors'] (no matching obsp graph)")
        del adata.uns["neighbors"]
        for key in ("distances", "connectivities"):
            if key in adata.obsp:
                del adata.obsp[key]

    t0 = time.time()

    work = adata
    if n_hvg and n_hvg > 0:
        # Select HVGs on a copy so the returned object keeps every gene — callers still need
        # the full matrix for the autocorrelation scores below.
        import numpy as np

        work = adata.copy()
        flavor = "seurat" if "log1p" in adata.uns else "seurat_v3"
        try:
            sc.pp.highly_variable_genes(work, n_top_genes=n_hvg, flavor=flavor)
        except Exception as exc:
            _log(f"  HVG selection with flavor={flavor} failed ({exc}); falling back to seurat")
            sc.pp.highly_variable_genes(work, n_top_genes=n_hvg, flavor="seurat")
        n_sel = int(work.var["highly_variable"].sum())
        work = work[:, work.var["highly_variable"]].copy()
        _log(f"  HVG selection ({flavor}): {n_sel} of {adata.n_vars} genes kept for the PCA")

    sc.pp.neighbors(work)
    sc.tl.umap(work)
    leiden_kwargs = dict(flavor="igraph", n_iterations=2, directed=False)
    if resolution is not None:
        leiden_kwargs["resolution"] = resolution
    sc.tl.leiden(work, **leiden_kwargs)

    if work is not adata:
        # Carry the labels, the embeddings, and the kNN graph back onto the full object.
        # The graph matters: without uns["neighbors"] + obsp, re-clustering the saved object
        # later fails with 'No "neighbors" in .uns' — the graph was built on the HVG subset
        # and would otherwise be discarded with it.
        adata.obs["leiden"] = work.obs["leiden"].values
        for key in ("X_pca", "X_umap"):
            if key in work.obsm:
                adata.obsm[key] = work.obsm[key]
        if "neighbors" in work.uns:
            adata.uns["neighbors"] = work.uns["neighbors"]
        for key in ("distances", "connectivities"):
            if key in work.obsp:
                adata.obsp[key] = work.obsp[key]

    n_clusters = adata.obs["leiden"].nunique()
    res_note = f", resolution={resolution}" if resolution is not None else " (default resolution)"
    _log(f"leiden computed in {_fmt(time.time() - t0)} -> {n_clusters} clusters{res_note}")
    sizes = adata.obs["leiden"].value_counts().tolist()
    _log(f"  cluster sizes: {sizes}")


def alias_cluster_scores_to_leiden(adata, cluster_key: str) -> None:
    """Copy `<cluster_key>_*` score entries to the `leiden_*` names the GeoJSON reads.

    Also mirrors the label column itself, because visium_to_geojson.py indexes the
    annotation blocks by `obs["leiden"]` cast to int.

    Note the labels necessarily become integers: the GeoJSON builder does
    `obs["leiden"].dropna().astype(int)`, so a categorical like "myocardium_ventricular"
    cannot survive the trip. The code -> name mapping is therefore also written to
    `uns["leiden_cluster_names"]` so the information is not lost, even though the current
    frontend renders the numbers ("Cluster 3") rather than the names.
    """
    for suffix in ("centrality_scores", "co_occurrence", "nhood_enrichment"):
        src, dst = f"{cluster_key}_{suffix}", f"leiden_{suffix}"
        if src in adata.uns:
            adata.uns[dst] = adata.uns[src]
            _log(f"aliased uns['{src}'] -> uns['{dst}']")
    if cluster_key != "leiden":
        cat = adata.obs[cluster_key].astype("category")
        codes = cat.cat.codes
        adata.obs["leiden"] = codes.astype(str).astype("category")
        mapping = {str(i): str(name) for i, name in enumerate(cat.cat.categories)}
        adata.uns["leiden_cluster_names"] = mapping
        _log(
            f"wrote obs['leiden'] from obs['{cluster_key}'] as integer codes "
            "(the GeoJSON casts leiden labels to int)"
        )
        for code, name in mapping.items():
            _log(f"    cluster {code} = {name}")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Backfill Squidpy scores into an already-processed AnnData.",
    )
    parser.add_argument("--adata", required=True, help="Path to the .h5ad to update.")
    parser.add_argument(
        "--out",
        default=None,
        help="Where to write. Defaults to updating --adata in place.",
    )
    parser.add_argument(
        "--cluster-key",
        default="leiden",
        help="obs column to group cluster scores by (default: leiden, computed if absent).",
    )
    parser.add_argument(
        "--alias-to-leiden",
        action="store_true",
        help="Also store the cluster scores under leiden_* so the GeoJSON picks them up.",
    )
    parser.add_argument("--interval", type=int, default=50, help="Co-occurrence distance bins.")
    parser.add_argument("--n-splits", type=int, default=None, help="Co-occurrence n_splits.")
    parser.add_argument(
        "--n-perms-nhood", type=int, default=1000, help="Neighbourhood-enrichment permutations."
    )
    parser.add_argument(
        "--n-perms-autocorr",
        type=int,
        default=None,
        help="Moran's I / Geary's C permutations. None uses the analytic approximation, "
        "which is what the pipeline does by default.",
    )
    parser.add_argument("--corr-method", default="fdr_bh", help="Multiple-testing correction.")
    parser.add_argument(
        "--skip",
        nargs="*",
        default=[],
        choices=["centrality", "co_occurrence", "nhood_enrichment", "moranI", "gearyC"],
        help="Scores to leave out.",
    )
    parser.add_argument(
        "--preserve-existing-leiden", dest="preserve_leiden", default=None, metavar="COLNAME",
        help="Before --alias-to-leiden overwrites obs['leiden'], copy the existing labels to "
             "COLNAME so an earlier unsupervised clustering is kept for comparison.",
    )
    parser.add_argument(
        "--n-hvg", dest="n_hvg", type=int, default=0,
        help="Restrict the clustering PCA to the N most variable genes (0 = all genes, the "
             "pipeline's own behaviour). ~2000 is the usual choice.",
    )
    parser.add_argument(
        "--resolution", type=float, default=None,
        help="Leiden resolution. Omit for scanpy's default (1.0). Lower gives fewer, larger "
             "clusters — try 0.3-0.5 when the tissue has only a handful of real states.",
    )
    parser.add_argument(
        "--force", action="store_true", help="Recompute even where results already exist."
    )
    parser.add_argument(
        "--backup",
        action="store_true",
        help="Copy the input to <name>.pre-backfill.h5ad before overwriting in place.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Report what is missing and would be computed, then exit without writing.",
    )
    args = parser.parse_args()

    if not os.path.exists(args.adata):
        raise SystemExit(f"not found: {args.adata}")

    out_path = args.out or args.adata
    in_place = os.path.abspath(out_path) == os.path.abspath(args.adata)

    _log(f"reading {args.adata}")
    t0 = time.time()
    adata = sc.read_h5ad(args.adata)
    _log(f"loaded {adata.n_obs} obs x {adata.n_vars} vars in {_fmt(time.time() - t0)}")

    wanted = {
        "centrality": f"{args.cluster_key}_centrality_scores",
        "co_occurrence": f"{args.cluster_key}_co_occurrence",
        "nhood_enrichment": f"{args.cluster_key}_nhood_enrichment",
        "moranI": "moranI",
        "gearyC": "gearyC",
    }
    todo = [
        name
        for name, uns_key in wanted.items()
        if name not in args.skip and (args.force or uns_key not in adata.uns)
    ]

    _log("current state:")
    for name, uns_key in wanted.items():
        state = "present" if uns_key in adata.uns else "MISSING"
        planned = "will compute" if name in todo else "skip"
        _log(f"    {name:18} uns['{uns_key}']  {state:8} -> {planned}")
    needs_leiden = args.cluster_key == "leiden" and "leiden" not in adata.obs.columns
    _log(f"    obs['{args.cluster_key}']{'':6} {'MISSING' if needs_leiden else 'present'}")

    if args.dry_run:
        _log("dry run, nothing written")
        return 0
    if not todo and not needs_leiden:
        _log("nothing to do (use --force to recompute)")
        return 0

    if args.cluster_key == "leiden":
        ensure_leiden(adata, force=args.force, n_hvg=args.n_hvg, resolution=args.resolution)
    elif args.cluster_key not in adata.obs.columns:
        raise SystemExit(
            f"obs['{args.cluster_key}'] not found. Available columns: "
            f"{sorted(adata.obs.columns)[:20]} ..."
        )

    cluster_scores = {"centrality", "co_occurrence", "nhood_enrichment"}
    if (cluster_scores & set(todo)) or {"moranI", "gearyC"} & set(todo):
        ensure_spatial_graph(adata, force=args.force)

    if "centrality" in todo:
        t0 = time.time()
        sq.gr.centrality_scores(adata, cluster_key=args.cluster_key, show_progress_bar=False)
        _log(f"centrality scores in {_fmt(time.time() - t0)}")

    if "co_occurrence" in todo:
        t0 = time.time()
        sq.gr.co_occurrence(
            adata,
            cluster_key=args.cluster_key,
            interval=args.interval,
            n_splits=args.n_splits,
            show_progress_bar=False,
        )
        _log(f"co-occurrence in {_fmt(time.time() - t0)}")

    if "nhood_enrichment" in todo:
        t0 = time.time()
        sq.gr.nhood_enrichment(
            adata,
            cluster_key=args.cluster_key,
            seed=42,
            n_perms=args.n_perms_nhood,
            show_progress_bar=False,
        )
        _log(f"neighbourhood enrichment in {_fmt(time.time() - t0)}")

    for mode, name in (("moran", "moranI"), ("geary", "gearyC")):
        if name not in todo:
            continue
        t0 = time.time()
        # transformation=True only when running the analytic variant, matching calc_scores.py
        sq.gr.spatial_autocorr(
            adata,
            mode=mode,
            seed=42,
            n_perms=args.n_perms_autocorr,
            transformation=args.n_perms_autocorr is None,
            corr_method=args.corr_method,
            show_progress_bar=False,
        )
        _log(f"{name} over {adata.n_vars} genes in {_fmt(time.time() - t0)}")

    if args.alias_to_leiden:
        if args.preserve_leiden and "leiden" in adata.obs.columns:
            if args.preserve_leiden in adata.obs.columns and not args.force:
                _log(f"obs['{args.preserve_leiden}'] already exists, not overwriting")
            else:
                adata.obs[args.preserve_leiden] = adata.obs["leiden"].values
                n = adata.obs[args.preserve_leiden].nunique()
                _log(f"preserved existing leiden ({n} clusters) as obs['{args.preserve_leiden}']")
        alias_cluster_scores_to_leiden(adata, args.cluster_key)

    if in_place and args.backup:
        backup = args.adata.replace(".h5ad", ".pre-backfill.h5ad")
        if os.path.exists(backup):
            _log(f"backup already exists, leaving it alone: {backup}")
        else:
            _log(f"backing up to {backup}")
            shutil.copy2(args.adata, backup)

    # AnnData rejects non-string obsm DataFrame column names on write; the pipeline does the
    # same coercion before saving (calc_scores.py's "ensure all obsm DataFrame column names
    # are strings" fix).
    for key, value in adata.obsm.items():
        if hasattr(value, "columns"):
            value.columns = value.columns.astype(str)

    _log(f"writing {out_path}")
    t0 = time.time()
    adata.write(out_path)
    _log(f"written in {_fmt(time.time() - t0)}")
    _log("now regenerate the GeoJSON so the frontend sees the new keys:")
    _log(
        f"    python3 backend/visium_to_geojson.py --adata {out_path} "
        f"--outpath <dataset_dir>/hexagons.geojson --data_type visium"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
