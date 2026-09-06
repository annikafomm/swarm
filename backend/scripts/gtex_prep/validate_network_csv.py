#!/usr/bin/env python3
"""
Check a GENIE3 or SPONGE network csv against SWARM's upload contract *before* uploading it.

Why this exists
---------------
Every way of getting these files wrong fails silently or nearly so:

  * wrong column names   -> SPONGE: "At least one of the sponge network files has an
                            unsupported file format." GENIE3: a `setnames` on absent columns.
  * wrong identifier type -> no error at all. The pipeline logs "With the chosen parameters no
                            regulons could be created." and produces an empty regulon set, so
                            the Regulatory Scores tab stays empty with no indication why.
  * SPONGE files swapped  -> tolerated by calc_scores.R (it detects and corrects the swap),
                            but worth knowing rather than guessing.

A processing run takes minutes to hours, so finding this out afterwards is expensive. This
does the checks in a second.

Contracts (from backend/calc_R_scores/calc_scores.R and utils.R):

  genie3      regulatoryGene, targetGene, weight            -- gene SYMBOLS
  sponge-int  geneA, geneB, df, cor, pcor, mscor, p.val, p.adj  -- ENSEMBL ids
  sponge-cent gene, degree, eigenvector, betweenness, page_rank  -- ENSEMBL ids

Examples
--------
    python3 validate_network_csv.py --genie3 genie3_heart_lv_network.csv
    python3 validate_network_csv.py --sponge-interactions sponge_heart_lv_interactionNetwork.csv \\
                                    --sponge-centralities sponge_heart_lv_networkAnalysis.csv
    # also confirm the ids actually intersect the dataset you will score
    python3 validate_network_csv.py --genie3 net.csv --against-adata backend/data/heart/plasmidpoop/adata_tg_scores.h5ad
"""

import argparse
import csv
import re
import sys

CONTRACTS = {
    "genie3": {
        "cols": ["regulatoryGene", "targetGene", "weight"],
        "id_cols": ["regulatoryGene", "targetGene"],
        "id_kind": "symbol",
        "upload_as": "--genie_network (Network Scores -> Upload Genie3 network)",
    },
    "sponge_interactions": {
        "cols": ["geneA", "geneB", "df", "cor", "pcor", "mscor", "p.val", "p.adj"],
        "id_cols": ["geneA", "geneB"],
        "id_kind": "ensembl",
        "upload_as": "--sponge_network (SPONGE 'networkinteractions')",
    },
    "sponge_centralities": {
        "cols": ["gene", "degree", "eigenvector", "betweenness", "page_rank"],
        "id_cols": ["gene"],
        "id_kind": "ensembl",
        "upload_as": "--sponge_analysis (SPONGE 'networkanalysis')",
    },
}

ENSEMBL_RE = re.compile(r"^ENSG\d{11}")


class Result:
    def __init__(self):
        self.errors, self.warnings = [], []

    def error(self, msg):
        self.errors.append(msg)

    def warn(self, msg):
        self.warnings.append(msg)


def sniff_id_kind(values) -> str:
    """Classify a sample of identifiers as ensembl / symbol / mixed."""
    ens = sum(1 for v in values if ENSEMBL_RE.match(v))
    if ens == 0:
        return "symbol"
    if ens >= len(values) * 0.9:
        return "ensembl"
    return "mixed"


def check(path: str, kind: str, res: Result, sample_ids: set) -> None:
    spec = CONTRACTS[kind]
    print(f"\n=== {kind}: {path} ===")
    try:
        fh = open(path, newline="", encoding="utf-8")
    except OSError as exc:
        res.error(f"{path}: cannot open ({exc})")
        return

    with fh:
        reader = csv.reader(fh)
        try:
            header = next(reader)
        except StopIteration:
            res.error(f"{path}: file is empty")
            return
        header = [h.strip() for h in header]
        print(f"  columns: {header}")

        missing = [c for c in spec["cols"] if c not in header]
        if missing:
            res.error(
                f"{path}: missing required column(s) {missing}. "
                f"Required: {spec['cols']}"
            )
            # Detect the classic swap so the message is actionable.
            other = {
                "sponge_interactions": "sponge_centralities",
                "sponge_centralities": "sponge_interactions",
            }.get(kind)
            if other and all(c in header for c in CONTRACTS[other]["cols"]):
                res.warn(
                    f"{path}: this looks like the {other} file instead. calc_scores.R does "
                    "detect and correct a swap, but pass them the right way round."
                )
            return
        extra = [c for c in header if c not in spec["cols"]]
        if extra:
            res.warn(f"{path}: extra columns present (harmless, ignored): {extra}")

        id_idx = [header.index(c) for c in spec["id_cols"]]
        ids, n_rows = [], 0
        for row in reader:
            n_rows += 1
            if len(row) < len(header):
                continue
            if len(ids) < 4000:
                ids.extend(row[i].strip() for i in id_idx)
        print(f"  data rows: {n_rows}")

    if n_rows == 0:
        res.error(f"{path}: header only, no data rows")
        return

    ids = [i for i in ids if i]
    if not ids:
        res.error(f"{path}: identifier columns are empty")
        return

    found = sniff_id_kind(ids)
    want = spec["id_kind"]
    print(f"  identifiers look like: {found} (e.g. {ids[0]!r}); contract wants: {want}")
    if found == "mixed":
        res.warn(f"{path}: mixed Ensembl/symbol identifiers — check the source matrix")
    elif found != want:
        res.error(
            f"{path}: identifiers are {found} but this network is matched against {want}. "
            + (
                "GENIE3 is matched against the AnnData var index (symbols)."
                if want == "symbol"
                else "SPONGE is matched against the AnnData ensembl column (ENSG ids)."
            )
            + " This fails SILENTLY — you get an empty regulon set, not an error."
        )

    if sample_ids:
        uniq = set(ids)
        overlap = len(uniq & sample_ids)
        pct = 100.0 * overlap / max(1, len(uniq))
        print(f"  overlap with dataset: {overlap} of {len(uniq)} network ids ({pct:.1f}%)")
        if overlap == 0:
            res.error(
                f"{path}: zero overlap with the dataset's genes — this network would produce "
                "nothing. Usually the wrong identifier type or the wrong species."
            )
        elif pct < 5:
            res.warn(f"{path}: only {pct:.1f}% of network ids appear in the dataset")


def load_dataset_ids(adata_path: str) -> set:
    """Collect both the var index and any Ensembl-looking var column from an .h5ad.

    Read with h5py rather than anndata so this stays fast on multi-hundred-MB files — we only
    need the var frame's string columns.
    """
    import h5py

    def _decode(values):
        out = []
        for v in values:
            if isinstance(v, bytes):
                out.append(v.decode("utf-8", "replace"))
            elif isinstance(v, str):
                out.append(v)
        return out

    ids = set()
    with h5py.File(adata_path, "r") as f:
        var = f["var"]
        idx_key = var.attrs.get("_index", "_index")
        if isinstance(idx_key, bytes):
            idx_key = idx_key.decode()
        for key in [idx_key] + [k for k in var.keys() if k != idx_key]:
            obj = var[key]
            try:
                if hasattr(obj, "keys"):
                    # Categorical column: anndata stores it as a group of categories + codes.
                    # Only the categories are needed for an overlap check, and skipping this
                    # branch is exactly why an annotated dataset can look unannotated —
                    # pandas turns a mostly-repeating string column into a categorical, so
                    # `ensembl_id` very often lands here rather than as a flat array.
                    if "categories" not in obj:
                        continue
                    decoded = _decode(obj["categories"][:200000])
                elif hasattr(obj, "shape") and len(obj.shape) == 1:
                    decoded = _decode(obj[:200000])
                else:
                    continue
            except Exception:
                continue
            if not decoded:
                continue
            # Keep the index plus anything that looks like Ensembl ids; skip numeric columns.
            sample = decoded[:200]
            if key == idx_key or sum(1 for d in sample if ENSEMBL_RE.match(d)) > 10:
                ids.update(d for d in decoded if d)
    return ids


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--genie3")
    ap.add_argument("--sponge-interactions", dest="sponge_int")
    ap.add_argument("--sponge-centralities", dest="sponge_cent")
    ap.add_argument(
        "--against-adata",
        default=None,
        help="Optional .h5ad to check identifier overlap against (catches species/id mismatches).",
    )
    args = ap.parse_args()

    if not any([args.genie3, args.sponge_int, args.sponge_cent]):
        ap.error("give at least one of --genie3 / --sponge-interactions / --sponge-centralities")

    sample_ids = set()
    if args.against_adata:
        print(f"reading var ids from {args.against_adata} ...")
        sample_ids = load_dataset_ids(args.against_adata)
        print(f"  {len(sample_ids)} identifiers available in the dataset")

    res = Result()
    if args.genie3:
        check(args.genie3, "genie3", res, sample_ids)
    if args.sponge_int:
        check(args.sponge_int, "sponge_interactions", res, sample_ids)
    if args.sponge_cent:
        check(args.sponge_cent, "sponge_centralities", res, sample_ids)

    if args.sponge_int and not args.sponge_cent:
        res.warn("SPONGE needs BOTH files; the upload form will not submit with only one")
    if args.sponge_cent and not args.sponge_int:
        res.warn("SPONGE needs BOTH files; the upload form will not submit with only one")

    print("\n" + "=" * 60)
    for w in res.warnings:
        print(f"WARN  {w}")
    for e in res.errors:
        print(f"ERROR {e}")
    if res.errors:
        print(f"\nFAILED with {len(res.errors)} error(s)")
        return 1
    print("\nPASSED" + (f" with {len(res.warnings)} warning(s)" if res.warnings else ""))
    for key, path in (("genie3", args.genie3), ("sponge_interactions", args.sponge_int),
                      ("sponge_centralities", args.sponge_cent)):
        if path:
            print(f"  {path}\n      -> upload as {CONTRACTS[key]['upload_as']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
