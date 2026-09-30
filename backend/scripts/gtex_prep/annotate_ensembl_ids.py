#!/usr/bin/env python3
"""
Add the `ensembl_id` and `feature_type` var columns that SPONGE scoring requires.

Why this exists
---------------
A SPONGE network is matched against the dataset by Ensembl id, not by symbol:

    ids <- var_df[, get(args$ensembl_col)]     # calc_scores.R
    rownames(expr_sp) <- ids

and calc_scores.R validates **both** of the columns the upload form exposes before it will run
anything:

    } else if (! args$ensembl_col %in% colnames(var_df)) {   # default "ensembl_id"
    } else if (! args$feature_col  %in% colnames(var_df)) {  # default "feature_type"

Datasets whose `var` carries only gene symbols therefore cannot use SPONGE at all — the run
aborts with "'ensembl_id' is not a column in adata.vars" no matter how good the network is.
The builtin `Heart (Multiome)` is exactly this case: its var columns are all Seurat/scanpy QC
metrics (`dispersions`, `vf_vst_counts_*`, ...) and the index holds symbols (`FAM87B`,
`LINC01409`).

`feature_type` matters beyond the existence check: SPONGE module building filters on it via
`--RNA_types` (default `lncRNA,protein_coding`), so a column of the wrong vocabulary silently
filters everything out.

Sources used, and why:
  * **symbol -> Ensembl** from the GTEx GCT's own `Name`/`Description` columns. Deliberately
    the same file the network was built from, so the two agree by construction — a separate
    annotation release could disagree on which Ensembl id a symbol maps to.
  * **Ensembl -> biotype** from `SPONGE::ensembl.df` (68,005 rows, `ensembl_gene_id` +
    `gene_biotype`), i.e. the vocabulary SPONGE itself filters against. Exported to tsv by a
    one-liner if not already present.

Nothing is overwritten: existing `ensembl_id` / `feature_type` columns are left alone unless
--force is given, and the input is only modified in place when you ask for it.

Examples
--------
    # inspect first
    python3 annotate_ensembl_ids.py --adata adata_tg_scores.h5ad \\
        --gct _downloads/gene_tpm_v10_heart_left_ventricle.gct.gz --dry-run

    # write in place, keeping a backup
    python3 annotate_ensembl_ids.py --adata adata_tg_scores.h5ad \\
        --gct _downloads/gene_tpm_v10_heart_left_ventricle.gct.gz \\
        --biotypes ensembl_biotypes.tsv --backup
"""

import argparse
import gzip
import os
import shutil
import subprocess
import sys


def _log(msg: str) -> None:
    print(f"[annotate] {msg}", flush=True)


def _open_text(path: str):
    return gzip.open(path, "rt", encoding="utf-8", errors="replace") if path.endswith(".gz") \
        else open(path, encoding="utf-8", errors="replace")


def symbol_to_ensembl(gct_path: str) -> dict:
    """Build {symbol: ENSG} from a GTEx GCT header block.

    Only the first two columns of each row are read, so this is cheap even on the full file.
    Where several Ensembl ids share a symbol the first is kept, matching how
    fetch_gtex_heart.py de-duplicates its symbol matrix — so the two stay consistent.
    """
    mapping = {}
    with _open_text(gct_path) as fh:
        fh.readline()  # "#1.2"
        fh.readline()  # dims
        header = fh.readline().rstrip("\n").split("\t")
        if header[:2] != ["Name", "Description"]:
            raise SystemExit(f"unexpected GCT header: {header[:4]}")
        for line in fh:
            tab1 = line.find("\t")
            if tab1 < 0:
                continue
            tab2 = line.find("\t", tab1 + 1)
            ens = line[:tab1].split(".")[0]
            sym = line[tab1 + 1:tab2 if tab2 > 0 else None].strip()
            if sym and sym not in mapping:
                mapping[sym] = ens
    return mapping


def ensembl_to_biotype(path: str | None) -> dict:
    """Build {ENSG: biotype}, exporting SPONGE::ensembl.df on demand if needed."""
    if path and os.path.exists(path):
        src = path
    else:
        src = path or "ensembl_biotypes.tsv"
        _log(f"exporting SPONGE::ensembl.df -> {src}")
        script = (
            'suppressPackageStartupMessages(library(SPONGE)); '
            'd <- SPONGE::ensembl.df[, c("ensembl_gene_id", "gene_biotype")]; '
            f'write.table(d, "{src}", sep="\\t", quote=FALSE, row.names=FALSE)'
        )
        res = subprocess.run(["Rscript", "-e", script], capture_output=True, text=True)
        if res.returncode != 0 or not os.path.exists(src):
            _log("could not export biotypes from SPONGE; feature_type will be left unset")
            _log(f"  {res.stderr.strip()[:200]}")
            return {}

    mapping = {}
    with open(src, encoding="utf-8", errors="replace") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            i_e, i_b = header.index("ensembl_gene_id"), header.index("gene_biotype")
        except ValueError:
            raise SystemExit(f"unexpected biotype columns in {src}: {header}")
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) > max(i_e, i_b):
                mapping[p[i_e].split(".")[0]] = p[i_b]
    return mapping


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--adata", required=True)
    ap.add_argument("--gct", required=True, help="GTEx GCT providing symbol -> Ensembl.")
    ap.add_argument("--biotypes", default=None,
                    help="tsv with ensembl_gene_id/gene_biotype; exported from SPONGE if absent.")
    ap.add_argument("--out", default=None, help="Defaults to updating --adata in place.")
    ap.add_argument("--ensembl-col", dest="ensembl_col", default="ensembl_id",
                    help="Column name to write [default %(default)s] — must match the upload "
                         "form's 'ensembl id (column name)' field.")
    ap.add_argument("--feature-col", dest="feature_col", default="feature_type",
                    help="Column name for biotype [default %(default)s] — must match the "
                         "form's 'feature (column name)' field.")
    ap.add_argument("--force", action="store_true", help="Overwrite existing columns.")
    ap.add_argument("--backup", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    if not os.path.exists(args.adata):
        raise SystemExit(f"not found: {args.adata}")

    import scanpy as sc

    _log(f"reading {args.adata}")
    adata = sc.read_h5ad(args.adata)
    _log(f"  {adata.n_obs} obs x {adata.n_vars} vars")
    _log(f"  existing var columns: {sorted(adata.var.columns)[:8]}"
         f"{' ...' if len(adata.var.columns) > 8 else ''}")

    for col in (args.ensembl_col, args.feature_col):
        if col in adata.var.columns and not args.force:
            _log(f"  var['{col}'] already present — use --force to overwrite")

    _log("building symbol -> Ensembl map from the GCT ...")
    sym2ens = symbol_to_ensembl(args.gct)
    _log(f"  {len(sym2ens)} symbols mapped")

    symbols = adata.var_names.astype(str)
    ens = [sym2ens.get(s, "") for s in symbols]
    hit = sum(1 for e in ens if e)
    _log(f"  matched {hit} of {len(symbols)} dataset genes ({100.0 * hit / len(symbols):.1f}%)")
    if hit == 0:
        raise SystemExit(
            "no dataset gene matched the GCT symbols — is the var index really gene symbols?"
        )

    ens2bio = ensembl_to_biotype(args.biotypes)
    if ens2bio:
        bio = [ens2bio.get(e, "") if e else "" for e in ens]
        bio_hit = sum(1 for b in bio if b)
        _log(f"  biotypes: {bio_hit} of {len(symbols)} genes annotated")
        from collections import Counter
        top = Counter(b for b in bio if b).most_common(4)
        _log(f"  most common: {', '.join(f'{k}={v}' for k, v in top)}")
    else:
        bio = None

    if args.dry_run:
        _log("dry run — would write:")
        _log(f"    var['{args.ensembl_col}'] ({hit} non-empty)")
        if bio:
            _log(f"    var['{args.feature_col}'] ({sum(1 for b in bio if b)} non-empty)")
        return 0

    if args.ensembl_col not in adata.var.columns or args.force:
        adata.var[args.ensembl_col] = ens
        _log(f"wrote var['{args.ensembl_col}']")
    if bio and (args.feature_col not in adata.var.columns or args.force):
        adata.var[args.feature_col] = bio
        _log(f"wrote var['{args.feature_col}']")

    out_path = args.out or args.adata
    if os.path.abspath(out_path) == os.path.abspath(args.adata) and args.backup:
        backup = args.adata.replace(".h5ad", ".pre-annotate.h5ad")
        if os.path.exists(backup):
            _log(f"backup already exists, leaving it: {backup}")
        else:
            _log(f"backing up to {backup}")
            shutil.copy2(args.adata, backup)

    # Same coercion the pipeline applies before writing: AnnData rejects non-string obsm
    # DataFrame column names.
    for key, value in adata.obsm.items():
        if hasattr(value, "columns"):
            value.columns = value.columns.astype(str)

    _log(f"writing {out_path}")
    adata.write(out_path)
    _log("done — SPONGE scoring can now resolve this dataset's genes")
    _log(f"  set the upload form's ensembl column to '{args.ensembl_col}'"
         f" and feature column to '{args.feature_col}'")
    return 0


if __name__ == "__main__":
    sys.exit(main())
