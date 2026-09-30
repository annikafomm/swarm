#!/usr/bin/env python3
"""
Write a human transcription-factor list for use as GENIE3's `regulators` argument.

Why this exists
---------------
GENIE3 fits one random forest per target gene, using the regulator set as predictors, so cost
scales as O(n_targets * n_regulators * n_samples). Handing it all ~58k GTEx genes as both
targets and regulators is not viable. Restricting regulators to known TFs is the standard
move and also makes the resulting network interpretable as regulator -> target.

The TF universe is taken from CollecTRI via decoupler — the same resource the LIANA+ stage
already falls back to (`dc.op.collectri()` in `calc_python_scores/calc_liana.py`), so the two
parts of SWARM agree on what counts as a TF. Requires network access on first call;
decoupler caches afterwards.

Output is one gene symbol per line, which is what `run_genie3.R --tf-list` expects. Symbols
(not Ensembl) because a GENIE3 network is matched against the AnnData var index — see the
note in run_genie3.R.

Examples
--------
    python3 make_tf_list.py --out tfs.txt
    # keep only TFs actually present in an expression matrix
    python3 make_tf_list.py --out tfs.txt --restrict-to expr_symbols_left_ventricle.csv
"""

import argparse
import csv
import sys


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--out", default="tfs.txt", help="Output file, one symbol per line.")
    ap.add_argument(
        "--restrict-to",
        default=None,
        help="Expression csv (first column = gene id). Only TFs present there are kept, "
             "so GENIE3 is never handed a regulator it has no data for.",
    )
    ap.add_argument(
        "--min-targets",
        type=int,
        default=0,
        help="Drop TFs with fewer than this many CollecTRI targets (0 = keep all).",
    )
    args = ap.parse_args()

    try:
        import decoupler as dc
    except ImportError:
        raise SystemExit("decoupler is not importable; it provides the CollecTRI resource")

    print("[tf-list] fetching CollecTRI ...", flush=True)
    net = dc.op.collectri()
    # CollecTRI columns are source/target/weight; older decoupler used 'source' too, so this
    # is stable, but fail loudly rather than silently producing an empty list.
    if "source" not in net.columns:
        raise SystemExit(f"unexpected CollecTRI columns: {list(net.columns)}")

    counts = net.groupby("source").size()
    if args.min_targets > 0:
        counts = counts[counts >= args.min_targets]
    tfs = sorted(counts.index.astype(str))
    print(f"[tf-list] CollecTRI: {len(net)} edges, {len(tfs)} TFs", flush=True)

    if args.restrict_to:
        present = set()
        with open(args.restrict_to, newline="", encoding="utf-8") as fh:
            reader = csv.reader(fh)
            next(reader, None)  # header
            for row in reader:
                if row:
                    present.add(row[0].strip())
        before = len(tfs)
        tfs = [t for t in tfs if t in present]
        print(
            f"[tf-list] restricted to {args.restrict_to}: {len(tfs)} of {before} TFs present",
            flush=True,
        )
        if not tfs:
            raise SystemExit(
                "no CollecTRI TF matched the expression matrix — check that its first column "
                "holds gene SYMBOLS (not Ensembl ids)"
            )

    with open(args.out, "w", encoding="utf-8") as fh:
        fh.write("\n".join(tfs) + "\n")
    print(f"[tf-list] wrote {args.out} ({len(tfs)} TFs)", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
