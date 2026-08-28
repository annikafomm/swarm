#!/usr/bin/env python3
"""
Make the heart measured-spatial object ready for regulatory scoring.

Why this is needed
------------------
The regulatory-score family reaches genes through two different namespaces, and the pipeline's
convention is that `var_names` holds SYMBOLS while Ensembl ids live in a var column:

    GENIE3 network   symbol-keyed   (regulatoryGene=CEBPA, targetGene=PLIN1)  -> matched on var_names
    SPONGE network   Ensembl-keyed  (geneA=ENSG00000125166)                    -> matched on --ensembl-col

The Tangram-projected object already follows that convention (symbol index + `ensembl_id` column),
which is how it satisfied both. The measured object is the other way round: Ensembl index with the
symbols in `var['SYMBOL']`. Handed to the scorer as-is, GENIE3 would match zero genes and the
regulons would come back empty rather than wrong -- a silent failure.

So this rewrites the measured object into the pipeline's convention:

    var_names        <- var['SYMBOL']
    var['ensembl_id'] <- the old Ensembl index (preserved, nothing is lost)
    var['feature_type'] <- copied from var['feature_types'] if present, for SPONGE's biotype filter

Duplicate symbols
-----------------
33,538 Ensembl ids do not map to 33,538 unique symbols. Rather than silently letting
`var_names_make_unique()` invent `GENE-1`/`GENE-2` names that match nothing in either network,
this keeps, for each symbol, the single Ensembl id with the highest total expression, and drops
the rest. That is deterministic, keeps the most informative copy, and is reported. Genes with an
empty/NA symbol are dropped for the same reason -- they cannot be addressed by a symbol network.

Gene filtering, and why scoring the measured object is slow either way
---------------------------------------------------------------------
`--min-cells` drops genes detected in fewer than N spots. It is worth having -- it roughly halves
the setup cost (regulon construction 1m17s -> 39s on this dataset) and 92.6% of GENIE3 targets and
617/634 regulators survive at the default of 10 -- but be clear about what it does NOT fix.

Scoring the measured object is intrinsically slow, and filtering genes barely helps:

    projected object (26,868 genes, dense)        AUCell 16 sec, VIPER 12 sec
    measured object  (33,538 genes, unfiltered)   killed at 5 h, no output
    measured object  (13,222 genes, filtered)     completed in ~6 h 15 m

The cause is the per-SPOT tie structure, not per-GENE sparsity, and the two are different axes.
Rank-based scorers (AUCell, GSVA, ssGSEA) rank all genes within each spot:

    projected spot: 26,868 genes, largest tie group        7
    measured  spot: 13,222 genes, largest tie group  ~11,600   (all the zeros)

A gene detected in 200 spots is still zero in the other 2,953, so dropping genes that are zero
*everywhere* leaves the per-spot tie groups nearly unchanged -- which is exactly what happened
here: 20,316 genes dropped, ~11,600 zeros per spot remaining, scoring still hours.

So this is a property of real sparse count data, and the pipeline's R scoring is effectively
tuned for dense projected input. Budget hours for the measured object, and do not expect a gene
filter to make it fast. It does complete, and the scores are worth it:
regulon Moran's I is +0.44 (AUCell) to +0.63 (GSVA) on the measured object versus ~0.00 on the
projected one -- see heart_projection_control_findings.md.

Usage
-----
    python3 prepare_measured_for_scoring.py \\
        --adata backend/data/heart/plasmidpoop/adata_spatial_measured_v2.h5ad \\
        --out   backend/data/heart/plasmidpoop/adata_st_scores_measured.h5ad
"""

import argparse
import os
import sys
import time


def log(msg):
    print(f"[prep-measured] {msg}", flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--adata", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--symbol-col", default="SYMBOL",
                    help="var column holding gene symbols [default SYMBOL].")
    ap.add_argument("--min-cells", type=int, default=10,
                    help="Drop genes detected in fewer than this many spots [default 10]. Set 0 "
                         "to keep everything, but read the docstring first: with 0 this dataset "
                         "did not finish scoring in 5 hours.")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    if not os.path.exists(args.adata):
        raise SystemExit(f"not found: {args.adata}")

    import numpy as np
    import scanpy as sc
    from scipy import sparse

    t0 = time.time()
    log(f"reading {args.adata}")
    adata = sc.read_h5ad(args.adata)
    log(f"    {adata.shape} in {time.time()-t0:.0f}s")

    if args.symbol_col not in adata.var.columns:
        raise SystemExit(
            f"var has no {args.symbol_col!r} column (available: {list(adata.var.columns)})"
        )

    ensembl = [str(x) for x in adata.var_names]
    symbols = [None if x is None else str(x) for x in adata.var[args.symbol_col]]
    blank = sum(1 for s in symbols if not s or s.lower() in ("nan", "none", "na"))
    log(f"symbols: {len(symbols)} total, {blank} blank/NA (will be dropped)")

    # per-gene total expression, to pick a winner among duplicate symbols
    X = adata.X
    totals = np.asarray(X.sum(axis=0)).ravel() if sparse.issparse(X) else \
        np.asarray(X).sum(axis=0)

    best = {}
    for i, s in enumerate(symbols):
        if not s or s.lower() in ("nan", "none", "na"):
            continue
        prev = best.get(s)
        if prev is None or totals[i] > totals[prev]:
            best[s] = i
    keep_idx = sorted(best.values())
    dropped_dup = len(symbols) - blank - len(keep_idx)
    log(f"unique symbols kept: {len(keep_idx)}   dropped as duplicate symbols: {dropped_dup}")
    log(f"    (for each symbol the Ensembl id with the highest total expression is kept)")

    if args.min_cells > 0:
        detected = np.asarray((X != 0).sum(axis=0)).ravel() if sparse.issparse(X) else \
            (np.asarray(X) != 0).sum(axis=0)
        before = len(keep_idx)
        keep_idx = [i for i in keep_idx if detected[i] >= args.min_cells]
        log(f"gene filter: detected in >= {args.min_cells} spots -> kept {len(keep_idx)} "
            f"of {before} ({before - len(keep_idx)} dropped)")
        log(f"    halves regulon-construction cost; does NOT make scoring fast -- the per-spot")
        log(f"    tie groups (~11,600 zeros per spot) dominate. See the module docstring.")

    if args.dry_run:
        log("dry run, nothing written")
        return 0

    out = adata[:, keep_idx].copy()
    old_ensembl = [ensembl[i] for i in keep_idx]
    new_symbols = [symbols[i] for i in keep_idx]

    out.var["ensembl_id"] = old_ensembl
    if "feature_types" in out.var.columns and "feature_type" not in out.var.columns:
        # SPONGE's biotype filter looks for `feature_type` (singular) by default.
        out.var["feature_type"] = out.var["feature_types"].astype(str).values
    out.var_names = new_symbols
    out.var.index.name = None

    # Record provenance so a later reader can tell this is a derived object.
    out.uns["prepared_for_scoring_from"] = os.path.basename(args.adata)
    out.uns["prepared_for_scoring_note"] = (
        "var_names reindexed from Ensembl to SYMBOL for symbol-keyed network matching; "
        "original Ensembl ids preserved in var['ensembl_id']. Duplicate symbols resolved by "
        "highest total expression."
    )

    # Same coercion the pipeline applies before writing.
    for key, value in out.obsm.items():
        if hasattr(value, "columns"):
            value.columns = value.columns.astype(str)

    log(f"writing {args.out}  ({out.shape})")
    out.write(args.out)
    log(f"done in {time.time()-t0:.0f}s")
    log("next: run scripts/backfill_regulatory_scores.py on this file with --description st")
    return 0


if __name__ == "__main__":
    sys.exit(main())
