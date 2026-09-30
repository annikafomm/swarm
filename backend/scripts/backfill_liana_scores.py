#!/usr/bin/env python3
"""
Backfill the LIANA+ score family into an already-processed dataset.

Why this exists
---------------
LIANA+ only runs during the original upload, gated on `scores.liana_plus`. A dataset
processed with that flag off (e.g. the builtin `Heart (Multiome)`) is missing all three
LIANA+ analyses and therefore the Ligand-Receptor, Cell Composition TF Activity, TF Activity
and Pathway Activity tabs.

**No user-supplied files are required.** LIANA+ is the one score family that can run purely
on defaults:

  * ligand-receptor  -> liana's built-in `consensus` resource, no file
  * TF activity      -> `dc.op.collectri()` when no GRN csv is given
  * pathway activity -> `dc.op.progeny()` when no pathway csv is given

The two `dc.op.*` calls download from OmniPath at runtime, so this needs outbound network
access. Check with `curl -I https://omnipathdb.org/about` before a long run.

This calls the same functions as the pipeline (`liana_folder.liana_bivariate`,
`liana_folder.liana_unsupervised`) with the same arguments as
`calc_python_scores/calc_liana.py:run_liana`, so a backfilled dataset is indistinguishable
from one processed with `scores.liana_plus` enabled. Splitting the stages out rather than
calling `run_liana` wholesale only buys `--skip` and per-stage timing.

Two behaviours inherited from the pipeline that are worth knowing before you run this:

  * **The ligand-receptor resource must match the organism.** LIANA's `consensus` resource
    is human-symbol based, `mouseconsensus` is mouse. This script picks one from
    `--organism` (default `auto`, inferred from symbol casing) and never modifies
    `adata.var_names`.

    Earlier versions of both this script and `run_liana` instead did
    `adata.var.index = adata.var.index.str.upper()` so mouse symbols would match the human
    resource. That mutation was permanent and destroyed the symbols
    (`Gapdh` -> `GAPDH`, `mt-Nd1` -> `MT-ND1`); it is not reversible by re-capitalising,
    since correct mouse casing follows no single rule. If you are working with a dataset
    processed before this was fixed, repair it with
    `scripts/restore_var_symbols.py` (needs an uncorrupted source file) before backfilling.
  * **`uns["liana_columns"]` is mandatory downstream.** The two cosine-similarity matrices
    are stored as bare arrays whose column names live only in that dict, and
    `visium_to_geojson.py` re-attaches them at build time. It is written here for whichever
    stages actually ran, so partial runs stay consistent.

After writing, regenerate the GeoJSON so the frontend sees the new keys:

    python3 backend/visium_to_geojson.py --adata <adata> \\
        --outpath <dataset_dir>/hexagons.geojson --data_type visium

Examples
--------
    python3 backend/scripts/backfill_liana_scores.py \\
        --adata backend/data/heart/plasmidpoop/adata_tg_scores.h5ad --dry-run

    python3 backend/scripts/backfill_liana_scores.py \\
        --adata backend/data/heart/plasmidpoop/adata_tg_scores.h5ad --backup
"""

import argparse
import os
import shutil
import sys
import time

# The pipeline modules import each other as top-level packages (`from liana_folder...`),
# so backend/ has to be importable regardless of where this is invoked from.
BACKEND_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for path in (BACKEND_DIR, os.path.join(BACKEND_DIR, "calc_python_scores")):
    if path not in sys.path:
        sys.path.insert(0, path)

import pandas as pd  # noqa: E402
import scanpy as sc  # noqa: E402


def _fmt(seconds: float) -> str:
    return f"{seconds:.1f}s" if seconds < 60 else f"{seconds / 60:.1f}min"


def _log(message: str) -> None:
    print(f"[liana-backfill] {message}", flush=True)


# Written by each stage; used to report what is missing and to verify afterwards.
STAGE_KEYS = {
    "ligand_receptor": {
        "obsm": [
            "ligand_receptor_cosine_similarity",
            "ligand_receptor_p_value",
            "ligand_receptor_category",
            "ligand_receptor_NMF_factors",
        ],
        "uns": ["ligand_receptor_global_scores"],
    },
    "cell_comp_tf_activity": {
        "obsm": [
            "cell_comp_tf_activity_cosine_similarity",
            "cell_comp_tf_activity_category",
            "tf_activity_score_ulm",
            "tf_activity_padj_ulm",
        ],
        "uns": ["cell_comp_tf_activity_global_scores"],
    },
    "pathway_activity": {
        "obsm": ["pathway_activity_score_mlm", "pathway_activity_padj_mlm"],
        "uns": [],
    },
}


def _stage_present(adata, stage: str) -> bool:
    keys = STAGE_KEYS[stage]
    return all(k in adata.obsm for k in keys["obsm"]) and all(
        k in adata.uns for k in keys["uns"]
    )


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Backfill LIANA+ scores into an already-processed AnnData.",
    )
    parser.add_argument("--adata", required=True, help="Path to the .h5ad to update.")
    parser.add_argument("--out", default=None, help="Defaults to updating --adata in place.")
    parser.add_argument(
        "--grn",
        default=None,
        help="TF-target GRN csv. Omit to download CollecTRI via decoupler.",
    )
    parser.add_argument(
        "--pathway-net",
        default=None,
        help="Pathway network csv. Omit to download PROGENy via decoupler.",
    )
    parser.add_argument(
        "--cell-comp-key",
        default="tangram_ct_pred",
        help="obsm key holding per-spot cell-type composition (default: tangram_ct_pred). "
        "The cell_comp_tf_activity stage is skipped when it is absent.",
    )
    parser.add_argument(
        "--skip",
        nargs="*",
        default=[],
        choices=sorted(STAGE_KEYS),
        help="Stages to leave out.",
    )
    parser.add_argument(
        "--organism",
        choices=("auto", "human", "mouse"),
        default="auto",
        help="Which ligand-receptor resource to use: human -> 'consensus', "
        "mouse -> 'mouseconsensus'. 'auto' infers it from gene-symbol casing.",
    )
    parser.add_argument("--force", action="store_true", help="Recompute existing results.")
    parser.add_argument(
        "--backup",
        action="store_true",
        help="Copy the input to <name>.pre-liana.h5ad before overwriting in place.",
    )
    parser.add_argument("--dry-run", action="store_true", help="Report and exit.")
    args = parser.parse_args()

    if not os.path.exists(args.adata):
        raise SystemExit(f"not found: {args.adata}")

    out_path = args.out or args.adata
    in_place = os.path.abspath(out_path) == os.path.abspath(args.adata)

    _log(f"reading {args.adata}")
    t0 = time.time()
    adata = sc.read_h5ad(args.adata)
    _log(f"loaded {adata.n_obs} obs x {adata.n_vars} vars in {_fmt(time.time() - t0)}")

    todo = []
    _log("current state:")
    for stage in sorted(STAGE_KEYS):
        present = _stage_present(adata, stage)
        blocked = ""
        if stage == "cell_comp_tf_activity" and args.cell_comp_key not in adata.obsm:
            blocked = f" (obsm['{args.cell_comp_key}'] missing)"
        planned = "skip"
        if stage not in args.skip and not blocked and (args.force or not present):
            planned = "will compute"
            todo.append(stage)
        _log(f"    {stage:24} {'present' if present else 'MISSING':8} -> {planned}{blocked}")

    if args.dry_run:
        _log("dry run, nothing written")
        return 0
    if not todo:
        _log("nothing to do (use --force to recompute)")
        return 0

    # Imported here so --dry-run and --help stay fast: these pull in liana, decoupler and
    # tensorflow, which together take ~20s and hit the network.
    from liana_folder.liana_bivariate import (  # noqa: E402
        cell_comp_tf_activity_similarity,
        ligand_receptor_relationships,
    )
    from liana_folder.liana_unsupervised import pathway_activities  # noqa: E402

    # Matches run_liana: the dict must exist before the bivariate stages append to it, and
    # must be preserved when it is already populated by an earlier partial run.
    adata.uns.setdefault("liana_columns", {})

    adata.var.index.name = None

    # Matches run_liana: pick the resource to fit the data, never rewrite the data to fit
    # the resource.
    organism = args.organism
    if organism == "auto":
        from calc_python_scores.calc_liana import looks_like_mouse  # noqa: E402

        organism = "mouse" if looks_like_mouse(adata) else "human"
        _log(f"organism auto-detected from symbol casing: {organism}")
    else:
        _log(f"organism specified: {organism}")
    resource_name = "mouseconsensus" if organism == "mouse" else "consensus"

    if "ligand_receptor" in todo:
        _log(f"ligand-receptor resource: {resource_name}")
        t0 = time.time()
        ligand_receptor_relationships(adata, resource_name=resource_name)
        _log(f"ligand-receptor relationships in {_fmt(time.time() - t0)}")

    if "cell_comp_tf_activity" in todo:
        import decoupler as dc

        t0 = time.time()
        grn = pd.read_csv(args.grn) if args.grn else dc.op.collectri(organism=organism)
        _log(f"GRN: {'file ' + args.grn if args.grn else f'CollecTRI ({organism})'} "
             f"({len(grn)} edges)")
        cell_comp_tf_activity_similarity(
            adata, grn, cell_comp_obsm_key=args.cell_comp_key
        )
        _log(f"cell-composition x TF activity in {_fmt(time.time() - t0)}")

    if "pathway_activity" in todo:
        import decoupler as dc

        t0 = time.time()
        pathway_net = (
            pd.read_csv(args.pathway_net)
            if args.pathway_net
            else dc.op.progeny(organism=organism)
        )
        _log(
            f"pathways: {'file ' + args.pathway_net if args.pathway_net else f'PROGENy ({organism})'} "
            f"({len(pathway_net)} rows)"
        )
        pathway_activities(adata, pathway_net)
        _log(f"pathway activity in {_fmt(time.time() - t0)}")

    _log("verifying written keys:")
    for stage in todo:
        for key in STAGE_KEYS[stage]["obsm"]:
            _log(f"    obsm['{key}']: {'OK' if key in adata.obsm else 'NOT WRITTEN'}")
        for key in STAGE_KEYS[stage]["uns"]:
            _log(f"    uns['{key}']:  {'OK' if key in adata.uns else 'NOT WRITTEN'}")
    _log(f"    uns['liana_columns'] keys: {sorted(adata.uns.get('liana_columns', {}))}")

    if in_place and args.backup:
        backup = args.adata.replace(".h5ad", ".pre-liana.h5ad")
        if os.path.exists(backup):
            _log(f"backup already exists, leaving it alone: {backup}")
        else:
            _log(f"backing up to {backup}")
            shutil.copy2(args.adata, backup)

    # AnnData refuses non-string obsm DataFrame column names on write; the pipeline applies
    # the same coercion before saving.
    for key, value in adata.obsm.items():
        if hasattr(value, "columns"):
            value.columns = value.columns.astype(str)

    _log(f"writing {out_path}")
    t0 = time.time()
    adata.write(out_path)
    _log(f"written in {_fmt(time.time() - t0)}")
    _log("now regenerate the GeoJSON:")
    _log(
        f"    python3 backend/visium_to_geojson.py --adata {out_path} "
        f"--outpath <dataset_dir>/hexagons.geojson --data_type visium"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
