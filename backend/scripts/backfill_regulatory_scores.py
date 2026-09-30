#!/usr/bin/env python3
"""
Backfill the regulatory-score family (AUCell / GSVA / ssGSEA / VIPER) onto an existing dataset,
given a GENIE3 and/or SPONGE network.

Why this exists
---------------
This is the one score family that cannot be recovered from the spatial data alone — it needs an
externally inferred regulatory network. A dataset uploaded without one has an empty
**Regulatory Scores** tab and no Global Regulatory Scores table, permanently. See
`backend/scripts/gtex_prep/` for producing the networks themselves.

Unlike the squidpy and LIANA+ backfills, this does **not** reimplement the scoring. It drives
the pipeline's own scripts in the same order `app.py` does, in a temporary work directory laid
out the way they expect:

    1. dump expr_info_<desc>/{expr.mtx,cells.txt,var.csv}   (as calc_scores.py does)
    2. Rscript calc_R_scores/calc_scores.R --dir <work> ...  -> Rscores_<desc>/*.csv|json
    3. python3 calc_python_scores/add_to_adata.py -indir <work> -Rscores
                                                            -> merged into obsm/uns
    4. copy the updated .h5ad back

so results are identical to a real pipeline run by construction. The scoring functions
(`calc_spongeffects_gsva`, VIPER regulon building, module creation) are all in
`calc_R_scores/utils.R` and are reached through the same entry point.

Prerequisites, both enforced by calc_scores.R before it does any work:

  * **A GENIE3 network** (`regulatoryGene,targetGene,weight`, keyed by gene **symbol**), and/or
    **both** SPONGE files (`geneA,geneB,df,cor,pcor,mscor,p.val,p.adj` +
    `gene,degree,eigenvector,betweenness,page_rank`, keyed by **Ensembl id**).
  * For SPONGE only: `var` must contain the ensembl and feature columns named by
    `--ensembl-col` / `--feature-col`. Most datasets do not — use
    `gtex_prep/annotate_ensembl_ids.py` to add them, and
    `gtex_prep/validate_network_csv.py --against-adata` to confirm the ids actually intersect
    before spending a run on it.

VIPER is GENIE3-only; there is no SPONGE variant.

Examples
--------
    python3 backfill_regulatory_scores.py --adata <adata> \\
        --genie3-network backend/networks/GENIE3/heart_gtex_lv/genie3_heart_lv_network.csv \\
        --dry-run

    python3 backfill_regulatory_scores.py --adata <adata> \\
        --genie3-network .../genie3_heart_lv_network.csv \\
        --sponge-network .../sponge_heart_lv_interactionNetwork.csv \\
        --sponge-analysis .../sponge_heart_lv_networkAnalysis.csv \\
        --backup
"""

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
import time

BACKEND_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# What calc_scores.R writes into Rscores_<desc>/, and therefore what should appear in obsm/uns.
EXPECTED_KEYS = {
    "genie3": [
        "aucell_scores_genie3",
        "spongeffects_GSVA_scores_genie3",
        "spongeffects_ssGSEA_scores_genie3",
        "viper_scores_genie3",
    ],
    "sponge": [
        "aucell_scores_sponge",
        "spongeffects_GSVA_scores_sponge",
        "spongeffects_ssGSEA_scores_sponge",
    ],
}


def _log(msg: str) -> None:
    print(f"[reg-backfill] {msg}", flush=True)


def _fmt(seconds: float) -> str:
    return f"{seconds:.1f}s" if seconds < 60 else f"{seconds / 60:.1f}min"


def dump_expr_info(adata, work_dir: str, description: str) -> str:
    """Write expr_info_<desc>/ exactly as calc_scores.py does.

    calc_scores.R reads the matrix from this directory, so the three filenames and the
    orientation (genes in rows after R's read, i.e. mmwrite of adata.X as-is) have to match.
    """
    import pandas as pd
    from scipy import io, sparse

    folder = os.path.join(work_dir, f"expr_info_{description}")
    os.makedirs(folder, exist_ok=True)
    io.mmwrite(os.path.join(folder, "expr.mtx"), sparse.csr_matrix(adata.X))
    pd.Series(adata.obs_names).to_csv(
        os.path.join(folder, "cells.txt"), index=False, header=False
    )
    adata.var.to_csv(os.path.join(folder, "var.csv"))
    return folder


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--adata", required=True)
    ap.add_argument("--out", default=None, help="Defaults to updating --adata in place.")
    ap.add_argument("--genie3-network", dest="genie3", default=None)
    ap.add_argument("--sponge-network", dest="sponge_net", default=None,
                    help="SPONGE interactionNetwork csv.")
    ap.add_argument("--sponge-analysis", dest="sponge_cent", default=None,
                    help="SPONGE networkAnalysis csv.")
    ap.add_argument("--description", choices=["st", "tg"], default=None,
                    help="Which pipeline slot this adata is. Inferred from the filename "
                         "(adata_tg_scores.h5ad -> tg) when omitted.")
    ap.add_argument("--algorithms", nargs="*", default=["aucell", "gsva", "ssgsea", "viper"],
                    choices=["aucell", "gsva", "ssgsea", "viper"],
                    help="Which scores to compute [default: all].")
    ap.add_argument("--ensembl-col", dest="ensembl_col", default="ensembl_id",
                    help="var column holding Ensembl ids (SPONGE only) [default %(default)s].")
    ap.add_argument("--feature-col", dest="feature_col", default="feature_type",
                    help="var column holding biotypes (SPONGE only) [default %(default)s].")
    ap.add_argument("--rna-types", dest="rna_types", default="lncRNA,protein_coding",
                    help="SPONGE biotype filter [default %(default)s].")
    ap.add_argument("--top-n", dest="top_n", type=int, default=100000,
                    help="GENIE3 edges to keep [default %(default)s].")
    ap.add_argument("--n-regulons", dest="n_regulons", type=int, default=20)
    ap.add_argument("--k-reg-genes", dest="k_reg_genes", type=int, default=20)
    ap.add_argument("--mscor", type=float, default=0.1)
    ap.add_argument("--padj", type=float, default=0.05)
    ap.add_argument("--max-modules", dest="max_modules", type=int, default=20)
    ap.add_argument("--n-cores", dest="n_cores", type=int, default=4)
    ap.add_argument("--keep-work-dir", dest="keep_work", action="store_true",
                    help="Leave the temp work dir in place for inspection.")
    ap.add_argument("--backup", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    if not os.path.exists(args.adata):
        raise SystemExit(f"not found: {args.adata}")

    has_genie = bool(args.genie3)
    has_sponge = bool(args.sponge_net and args.sponge_cent)
    if args.sponge_net and not args.sponge_cent or args.sponge_cent and not args.sponge_net:
        raise SystemExit(
            "SPONGE needs BOTH --sponge-network and --sponge-analysis "
            "(calc_scores.R validates the pair)"
        )
    if not has_genie and not has_sponge:
        raise SystemExit(
            "give --genie3-network and/or both SPONGE files — there is nothing to score against"
        )
    for path in filter(None, (args.genie3, args.sponge_net, args.sponge_cent)):
        if not os.path.exists(path):
            raise SystemExit(f"not found: {path}")

    description = args.description
    if description is None:
        description = "tg" if os.path.basename(args.adata).endswith("tg_scores.h5ad") else "st"
    _log(f"pipeline slot: {description} "
         f"({'Tangram-projected' if description == 'tg' else 'spatial'})")

    import scanpy as sc

    _log(f"reading {args.adata}")
    t0 = time.time()
    adata = sc.read_h5ad(args.adata)
    _log(f"  {adata.n_obs} obs x {adata.n_vars} vars in {_fmt(time.time() - t0)}")

    # Report what is already there, and what this run would add.
    planned = []
    if has_genie:
        planned += [k for k in EXPECTED_KEYS["genie3"]
                    if any(a in k.lower() or a == "viper" and "viper" in k for a in args.algorithms)]
    if has_sponge:
        planned += EXPECTED_KEYS["sponge"]
    planned = [k for k in dict.fromkeys(planned)
               if any(alg in k.lower() for alg in args.algorithms)]
    _log("current state:")
    for key in planned:
        _log(f"    obsm['{key}']: {'present' if key in adata.obsm else 'MISSING'}")

    # SPONGE's hard prerequisites, checked before doing minutes of work.
    if has_sponge:
        missing_cols = [c for c in (args.ensembl_col, args.feature_col)
                        if c not in adata.var.columns]
        if missing_cols:
            raise SystemExit(
                f"SPONGE needs var column(s) {missing_cols}, which this dataset lacks. "
                "calc_scores.R aborts on that check. Add them with "
                "gtex_prep/annotate_ensembl_ids.py, or drop the SPONGE arguments."
            )
        n_ens = int((adata.var[args.ensembl_col].astype(str).str.len() > 0).sum())
        _log(f"  var['{args.ensembl_col}']: {n_ens} of {adata.n_vars} non-empty")
        if n_ens == 0:
            raise SystemExit(f"var['{args.ensembl_col}'] is entirely empty")

    if args.dry_run:
        _log("dry run — no work dir created, nothing run")
        return 0

    work_dir = tempfile.mkdtemp(prefix="reg_backfill_")
    _log(f"work dir: {work_dir}")
    try:
        # The R script and add_to_adata both key off the *_scores.h5ad filename convention.
        staged_name = f"adata_{description}_scores.h5ad"
        staged = os.path.join(work_dir, staged_name)
        _log(f"staging as {staged_name} (the name both pipeline steps key off)")
        t0 = time.time()
        adata.write(staged)
        _log(f"  staged in {_fmt(time.time() - t0)}")

        t0 = time.time()
        _log("dumping expression matrix for R ...")
        dump_expr_info(adata, work_dir, description)
        _log(f"  done in {_fmt(time.time() - t0)}")
        del adata  # free ~GBs before the subprocesses run

        log_file = os.path.join(work_dir, "backfill.log")
        r_cmd = [
            "Rscript", os.path.join(BACKEND_DIR, "calc_R_scores", "calc_scores.R"),
            "--dir", work_dir,
            "--log", log_file,
            "--ensembl_col", args.ensembl_col,
            "--feature_col", args.feature_col,
            "--n_cores", str(args.n_cores),
        ]
        # calc_scores.R picks which expr_info_* slot to read from this flag, not from any
        # path argument: `--tangram` -> expr_info_tg, absent -> expr_info_st (see its
        # dir_expr assignment and the two mutually-exclusive existence checks).
        if description == "tg":
            r_cmd.append("--tangram")

        for alg in args.algorithms:
            r_cmd.append(f"--{alg}")
        if has_genie:
            r_cmd += [
                "--genie_network", os.path.abspath(args.genie3),
                "--top_n", str(args.top_n),
                "--n_regulons", str(args.n_regulons),
                "--k_reg_genes", str(args.k_reg_genes),
            ]
        if has_sponge:
            r_cmd += [
                "--sponge_network", os.path.abspath(args.sponge_net),
                "--sponge_analysis", os.path.abspath(args.sponge_cent),
                "--mscor", str(args.mscor),
                "--padj", str(args.padj),
                "--max_modules", str(args.max_modules),
                "--RNA_types", args.rna_types,
            ]

        _log("running calc_scores.R (the long step) ...")
        _log(f"  {' '.join(r_cmd[:3])} ... {' '.join(a for a in r_cmd if a.startswith('--'))}")
        t0 = time.time()
        # cwd matters: the pipeline invokes these from frontend/, and calc_scores.R sources
        # utils.R by a path relative to its own location, so run from backend/ to be safe.
        res = subprocess.run(r_cmd, cwd=BACKEND_DIR, capture_output=True, text=True)
        _log(f"  calc_scores.R exited {res.returncode} after {_fmt(time.time() - t0)}")
        tail = (res.stdout or "").strip().splitlines()[-12:]
        for line in tail:
            _log(f"  R| {line}")
        if res.returncode != 0:
            for line in (res.stderr or "").strip().splitlines()[-15:]:
                _log(f"  R!| {line}")
            raise SystemExit("calc_scores.R failed; nothing written")

        scores_dir = os.path.join(work_dir, f"Rscores_{description}")
        if not os.path.isdir(scores_dir) or not os.listdir(scores_dir):
            # calc_scores.R logs rather than exits non-zero for several soft failures
            # ("no regulons could be created", bad ensembl column, ...), so check the output.
            if os.path.exists(log_file):
                for line in open(log_file, encoding="utf-8", errors="replace").read().splitlines()[-15:]:
                    _log(f"  log| {line}")
            raise SystemExit(
                f"no score files in {scores_dir}. calc_scores.R often reports this kind of "
                "problem in its log rather than by failing — the usual cause is a network "
                "whose identifiers do not match the dataset (check with "
                "gtex_prep/validate_network_csv.py --against-adata)."
            )
        produced = sorted(os.listdir(scores_dir))
        _log(f"R produced {len(produced)} file(s): {', '.join(produced)}")

        _log("merging scores back into the AnnData ...")
        t0 = time.time()
        merge_cmd = [
            "python3", os.path.join(BACKEND_DIR, "calc_python_scores", "add_to_adata.py"),
            "-indir", work_dir,
            "-log", log_file,
            "-Rscores",
        ]
        res = subprocess.run(merge_cmd, cwd=BACKEND_DIR, capture_output=True, text=True)
        _log(f"  add_to_adata exited {res.returncode} after {_fmt(time.time() - t0)}")
        for line in (res.stdout or "").strip().splitlines()[-8:]:
            _log(f"  py| {line}")
        if res.returncode != 0:
            for line in (res.stderr or "").strip().splitlines()[-15:]:
                _log(f"  py!| {line}")
            raise SystemExit("add_to_adata.py failed; nothing written")

        # Verify before overwriting the user's file.
        merged = sc.read_h5ad(staged)
        got = [k for k in planned if k in merged.obsm]
        _log(f"verified {len(got)} of {len(planned)} expected obsm key(s) present:")
        for key in planned:
            _log(f"    obsm['{key}']: {'OK' if key in merged.obsm else 'not written'}")
        gene_sets = [k for k in ("genie_genesets", "sponge_genesets") if k in merged.uns]
        if gene_sets:
            _log(f"  uns gene sets: {', '.join(f'{k} ({len(merged.uns[k])})' for k in gene_sets)}")
        if not got:
            raise SystemExit("no expected keys were written; leaving the input untouched")
        del merged

        out_path = args.out or args.adata
        if os.path.abspath(out_path) == os.path.abspath(args.adata) and args.backup:
            backup = args.adata.replace(".h5ad", ".pre-regulatory.h5ad")
            if os.path.exists(backup):
                _log(f"backup already exists, leaving it: {backup}")
            else:
                _log(f"backing up to {backup}")
                shutil.copy2(args.adata, backup)

        _log(f"writing {out_path}")
        shutil.move(staged, out_path)
        _log("done. now regenerate the GeoJSON so the frontend sees the new scores:")
        _log(f"    python3 backend/visium_to_geojson.py --adata {out_path} "
             f"--outpath <dataset_dir>/hexagons.geojson --data_type visium")

        # The pipeline also drops filtered copies of the networks next to the adata; keep them
        # since the frontend's network views read them via the dataset's *_network_path fields.
        for name in (f"genie_network_filtered_{description}.csv",
                     f"sponge_network_filtered_{description}.csv"):
            src = os.path.join(work_dir, name)
            if os.path.exists(src):
                dst = os.path.join(os.path.dirname(os.path.abspath(out_path)), name)
                shutil.copy2(src, dst)
                _log(f"  also wrote {dst}")
    finally:
        if args.keep_work:
            _log(f"work dir kept: {work_dir}")
        else:
            shutil.rmtree(work_dir, ignore_errors=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
