#!/usr/bin/env python3
"""
Restore original gene symbols to an AnnData whose `var` index was destructively upper-cased.

Why this exists
---------------
`calc_python_scores/calc_liana.py` upper-cases `adata.var.index` in place before running
LIANA+, because the ligand-receptor stage historically requested the human-symbol `consensus`
resource. The mutation is permanent and the caller writes the same object back to disk, so any
MOUSE dataset that went through LIANA+ has had its symbols rewritten: `Xkr4` -> `XKR4`,
`Gm1992` -> `GM1992`, `mt-Nd1` -> `MT-ND1`.

That breaks every join against a mouse annotation (MGI, Ensembl mouse, mouse GO sets,
mouse TargetScan) and makes the symbols look human, so a join can silently match the wrong
species' gene instead of failing loudly.

**It cannot be repaired by re-capitalising.** `str.upper()` is not invertible, because correct
mouse casing is not one pattern. Neither `.title()` nor `.capitalize()` recovers it, and they
fail on different subsets:

    uppercased        .title()        .capitalize()   correct
    GM1992            Gm1992  ok      Gm1992  ok      Gm1992
    4732440D04RIK     ...Rik  ok      ...rik  WRONG   4732440D04Rik
    RB1CC1            Rb1Cc1  WRONG   Rb1cc1  ok      Rb1cc1
    PARD3BOS1         Pard3Bos1 WRONG Pard3bos1 ok    Pard3bos1
    MT-ND1            Mt-Nd1  WRONG   Mt-nd1  WRONG   mt-Nd1
    AI597479          Ai597479 WRONG  Ai597479 WRONG  AI597479

and the uppercased string carries no information about which rule applies. On the visual-cortex
builtin ~14,700 of 31,053 symbols fall in at-risk classes (7,518 `Gm#####`, 4,596 with mid-name
digits, 2,174 Riken `...Rik`, 123 accession-style, 13 mitochondrial `mt-`).

So this script does not guess. It copies the symbols positionally from an uncorrupted source
file, and refuses to write unless the two objects provably correspond:

  1. same number of vars
  2. elementwise `[s.upper() for s in source] == corrupted`  (i.e. identical order)
  3. no uppercase collisions within the source (which would make any mapping ambiguous)

If check 2 fails it falls back to a case-insensitive *mapping* (still exact, just order-
independent), and refuses entirely if that is ambiguous.

Derived tables
--------------
`var_names` is not the only place the mangled symbols landed. Score stages that ran after
LIANA wrote per-gene result tables into `uns`, keyed by the symbols as they stood at the
time, so those indices are upper-cased too:

    uns/moranI, uns/gearyC                          one row per gene (31,053 for visual_cortex)
    uns/{viper,aucell,spongeffects_*}_scores_*      one row per regulon/TF

Those are repaired in the same pass, by mapping each index entry through
`{v.upper(): v for v in restored_var_names}`. That is exact rather than inferred, and only
runs when the restored var index has no uppercase collisions (already required above), so a
symbol can only ever map back to one original. Pass `--skip-uns` to leave them alone.

Not repaired: `uns/ligand_receptor_global_scores` and `uns/liana_columns["ligand_receptor"]`.
Those hold `LIGAND^RECEPTOR` pair names from LIANA's **human** `consensus` resource, so they
are not mangled mouse symbols and there is nothing to map them back to. To get mouse-native
pairs, re-run the ligand-receptor stage against `mouseconsensus`:

    python3 backend/scripts/backfill_liana_scores.py --adata <file> \
        --organism mouse --skip cell_comp_tf_activity pathway_activity --force

Examples
--------
    python3 restore_var_symbols.py \\
        --adata backend/data/visual_cortex/dumpase1/adata_st_scores.h5ad \\
        --source backend/data/visual_cortex/visual_cortex_visium_1_prepro_final.h5ad --dry-run

    python3 restore_var_symbols.py --adata <target> --source <original> --backup
"""

import argparse
import os
import shutil
import sys
from collections import Counter


def _log(msg: str) -> None:
    print(f"[restore-symbols] {msg}", flush=True)


def var_names(path: str, column: str | None = None):
    """Read a var index -- or a var column -- without loading the matrix.

    `column` exists because not every source keys its var by symbol. The heart builtin's
    spatial input is indexed by Ensembl ID (`ENSG00000243485`) and carries the symbols in
    `var["SYMBOL"]`, which is the authoritative HGNC casing for that dataset.
    """
    import h5py

    with h5py.File(path, "r") as f:
        v = f["var"]
        if column is None:
            key = v.attrs.get("_index", "_index")
            if isinstance(key, bytes):
                key = key.decode()
            node = v[key]
        else:
            if column not in v:
                raise SystemExit(
                    f"{path}: var has no column {column!r} (available: {sorted(v)})"
                )
            node = v[column]

        # Categoricals are stored as a group of categories + integer codes.
        if isinstance(node, h5py.Group):
            cats = [
                x.decode() if isinstance(x, bytes) else str(x)
                for x in node["categories"][:]
            ]
            return [cats[c] if c >= 0 else "" for c in node["codes"][:]]
        return [
            x.decode() if isinstance(x, bytes) else str(x)
            for x in node[:]
        ]


def uns_table_indices(path: str):
    """Return {uns_key: [index values]} for every DataFrame stored under uns.

    Read straight from the h5 so --dry-run stays cheap. AnnData writes a DataFrame as a
    group carrying an `_index` attribute naming the column that holds the index.
    """
    import h5py

    out = {}
    with h5py.File(path, "r") as f:
        if "uns" not in f:
            return out
        for key in f["uns"]:
            node = f["uns"][key]
            if not isinstance(node, h5py.Group):
                continue
            idx_key = node.attrs.get("_index")
            if idx_key is None:
                continue
            if isinstance(idx_key, bytes):
                idx_key = idx_key.decode()
            try:
                out[key] = [
                    x.decode() if isinstance(x, bytes) else str(x)
                    for x in node[idx_key][:]
                ]
            except Exception:
                # Non-string index (e.g. integer-keyed) — nothing to restore.
                continue
    return out


# uns entries that are {regulon_name: [target gene symbols]} rather than DataFrames.
# Written by calc_scores.R (`genie_genesets.json`) and loaded by add_to_adata.py; both the
# keys and the values are gene symbols, so both were mangled.
GENESET_KEYS = ("genie_genesets", "sponge_genesets")


def uns_geneset_members(path: str):
    """Return {uns_key: {regulon: [targets]}} for the geneset dicts stored under uns."""
    import h5py

    out = {}
    with h5py.File(path, "r") as f:
        if "uns" not in f:
            return out
        for key in GENESET_KEYS:
            node = f["uns"].get(key)
            if not isinstance(node, h5py.Group):
                continue
            members = {}
            for reg in node:
                try:
                    members[reg] = [
                        x.decode() if isinstance(x, bytes) else str(x)
                        for x in node[reg][:]
                    ]
                except Exception:
                    members[reg] = None
            out[key] = members
    return out


def plan_geneset_restore(path: str, restored: list):
    """Which geneset dicts can be repaired, and how many names change.

    Returns {uns_key: {"regulons": {old: new}, "targets": {old: new}, "n_targets": int}}.
    Regulon keys and target symbols are mapped through the same collision-free
    upper -> original table as everything else, so this is exact, not inferred.
    """
    upper_to_orig = {v.upper(): v for v in restored}
    plan = {}
    for key, members in uns_geneset_members(path).items():
        regulons, targets, n_targets, unresolved = {}, {}, 0, 0
        for reg, tgts in members.items():
            orig_reg = upper_to_orig.get(reg.upper())
            if orig_reg is None:
                unresolved += 1
            elif orig_reg != reg:
                regulons[reg] = orig_reg
            for t in tgts or []:
                n_targets += 1
                orig_t = upper_to_orig.get(t.upper())
                if orig_t is None:
                    unresolved += 1
                elif orig_t != t:
                    targets[t] = orig_t
        if regulons or targets:
            plan[key] = {
                "regulons": regulons,
                "targets": targets,
                "n_targets": n_targets,
                "unresolved": unresolved,
            }
    return plan


# obsm matrices whose COLUMNS are regulon/TF names, so they were mangled too. Matched by the
# suffix the pipeline uses (`visium_to_geojson.py` identifies them the same way) rather than by
# sniffing the contents, because several other obsm tables have columns that look like gene
# symbols but are not:
#   pathway_activity_*   PROGENy pathway names — one is literally "EGFR"
#   proportions_class    cell-type labels — Sst, Vip, Pvalb, Lamp5, Sncg, Meis2 are all real
#                        mouse genes as well as class names
# Rewriting either of those would silently relabel scores with the wrong meaning, which is
# worse than leaving the casing wrong. The all-or-nothing check below is the second guard.
OBSM_REGULON_SUFFIXES = ("_genie3", "_sponge")


def obsm_columns(path: str):
    """Return {obsm_key: [column names]} for DataFrame-valued obsm entries."""
    import h5py

    out = {}
    with h5py.File(path, "r") as f:
        if "obsm" not in f:
            return out
        for key in f["obsm"]:
            node = f["obsm"][key]
            if not isinstance(node, h5py.Group):
                continue
            cols = node.attrs.get("column-order")
            if cols is None:
                continue
            out[key] = [
                c.decode() if isinstance(c, bytes) else str(c) for c in cols
            ]
    return out


def plan_obsm_restore(path: str, restored: list):
    """Which regulon-score obsm tables can be repaired. Returns (plan, skipped)."""
    upper_to_orig = {v.upper(): v for v in restored}
    plan, skipped = {}, []
    for key, cols in obsm_columns(path).items():
        if not key.endswith(OBSM_REGULON_SUFFIXES):
            continue
        hits = sum(1 for c in cols if c.upper() in upper_to_orig)
        if hits < len(cols):
            skipped.append((key, len(cols), f"only {hits}/{len(cols)} resolvable"))
            continue
        mapping = {
            c: upper_to_orig[c.upper()] for c in cols if upper_to_orig[c.upper()] != c
        }
        if mapping:
            plan[key] = mapping
    return plan, skipped


def plan_uns_restore(path: str, restored: list):
    """Which uns tables can be repaired from the restored symbols, and how many rows change.

    Returns (plan, unresolvable) where plan maps uns_key -> {old: new} and `unresolvable`
    lists tables whose index is not gene symbols at all (LR pair names, cluster ids), which
    are reported and left untouched rather than treated as an error.
    """
    upper_to_orig = {v.upper(): v for v in restored}
    plan, unresolvable = {}, []
    for key, idx in uns_table_indices(path).items():
        if not idx:
            continue
        hits = sum(1 for x in idx if x.upper() in upper_to_orig)
        if hits == 0:
            unresolvable.append((key, len(idx), idx[:2]))
            continue
        if hits < len(idx):
            # A partially-matching index would leave the table half-repaired, which is
            # worse than not touching it — bail loudly instead of guessing.
            unresolvable.append((key, len(idx), [f"only {hits}/{len(idx)} resolvable"]))
            continue
        mapping = {x: upper_to_orig[x.upper()] for x in idx if upper_to_orig[x.upper()] != x}
        if mapping:
            plan[key] = mapping
    return plan, unresolvable


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--adata", required=True, help="File whose var index is upper-cased.")
    ap.add_argument("--source", required=True,
                    help="Uncorrupted file to take the original symbols from.")
    ap.add_argument("--out", default=None, help="Defaults to updating --adata in place.")
    ap.add_argument("--backup", action="store_true",
                    help="Copy the input to <name>.pre-symbol-restore.h5ad first.")
    ap.add_argument("--source-column", default=None,
                    help="Take the original symbols from this var column of --source instead "
                         "of its index (e.g. SYMBOL for the heart spatial input, which is "
                         "indexed by Ensembl ID).")
    ap.add_argument("--allow-unmatched", action="store_true",
                    help="Permit target symbols that are absent from --source, leaving them "
                         "unchanged instead of refusing. Use when the source is a different "
                         "annotation build rather than the same gene set -- the count is "
                         "always reported so you can judge whether it looks like the wrong "
                         "file.")
    ap.add_argument("--skip-uns", action="store_true",
                    help="Do not repair gene-keyed tables in uns (moranI, gearyC, regulon "
                         "score tables). See the module docstring.")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    for p in (args.adata, args.source):
        if not os.path.exists(p):
            raise SystemExit(f"not found: {p}")

    target = var_names(args.adata)
    source = [g for g in var_names(args.source, args.source_column) if g]
    _log(f"target  {os.path.basename(args.adata)}: {len(target)} vars, e.g. {target[:4]}")
    src_desc = f"{os.path.basename(args.source)}"
    if args.source_column:
        src_desc += f" [var['{args.source_column}']]"
    _log(f"source  {src_desc}: {len(source)} symbols, e.g. {source[:4]}")

    already_restored = not all(s == s.upper() for s in target[:200])
    if already_restored:
        # Either this file was never mangled, or var was repaired by an earlier run while the
        # derived uns tables still carry the upper-cased symbols. Both are fine to continue
        # from: the restore below is a no-op for var and the uns pass does the remaining work.
        # Guard against the wrong --source by requiring the two indices to agree case-
        # insensitively, which is the same correspondence the mangled case relies on.
        if {t.upper() for t in target} != {v.upper() for v in source}:
            _log("target is NOT upper-cased and its symbols do not match --source "
                 "case-insensitively — refusing, this looks like the wrong file pair")
            return 1
        _log("var index is already in original casing (nothing mangled, or a previous run "
             "repaired it) — continuing to check the derived uns tables")

    collisions = {k: n for k, n in Counter(s.upper() for s in source).items() if n > 1}
    if collisions:
        _log(f"REFUSING: {len(collisions)} uppercase collisions in the source, e.g. "
             f"{list(collisions)[:5]} — the mapping would be ambiguous")
        return 1
    _log("no uppercase collisions in the source (mapping is unambiguous)")

    strategy = None
    if already_restored and source == target:
        strategy = "already-restored"
        restored = list(target)
        _log("var needs no change -> uns-only pass")
    elif len(source) == len(target) and [s.upper() for s in source] == target:
        strategy = "positional"
        restored = list(source)
        _log("order matches elementwise -> positional restore (exact, no inference)")
    else:
        # Order differs (e.g. the pipeline subset or reordered genes). A case-insensitive
        # lookup is still exact because collisions were just ruled out.
        upper_to_orig = {s.upper(): s for s in source}
        missing = [t for t in target if t.upper() not in upper_to_orig]
        if missing and not args.allow_unmatched:
            _log(f"REFUSING: {len(missing)} of {len(target)} target symbols are absent from "
                 f"the source, e.g. {missing[:5]}")
            _log("    If the source is the same gene set, this means it is the wrong file. If "
                 "it is a different annotation build, pass --allow-unmatched to leave those "
                 "symbols as they are.")
            return 1
        strategy = "mapping"
        restored = [upper_to_orig.get(t.upper(), t) for t in target]
        _log("order differs -> case-insensitive mapping restore (exact, collision-free)")
        if missing:
            _log(f"    {len(missing)} of {len(target)} symbols absent from the source, left "
                 f"unchanged (--allow-unmatched), e.g. {missing[:4]}")

    changed = sum(1 for a, b in zip(target, restored) if a != b)
    _log(f"would change {changed} of {len(target)} symbols in var ({strategy})")
    for a, b in list(zip(target, restored))[:6]:
        if a != b:
            _log(f"    {a}  ->  {b}")

    uns_plan, uns_skipped, geneset_plan, obsm_plan = {}, [], {}, {}
    if not args.skip_uns:
        uns_plan, uns_skipped = plan_uns_restore(args.adata, restored)
        _log("gene-keyed tables in uns:")
        for key, mapping in sorted(uns_plan.items()):
            sample = ", ".join(f"{k}->{v}" for k, v in list(mapping.items())[:2])
            _log(f"    {key:44} {len(mapping):6} rows  ({sample})")
        for key, n, why in uns_skipped:
            _log(f"    {key:44} {n:6} rows  LEFT AS-IS ({', '.join(map(str, why))})")
        if not uns_plan:
            _log("    (none need repair)")

        geneset_plan = plan_geneset_restore(args.adata, restored)
        if geneset_plan:
            _log("regulon geneset dicts in uns:")
            for key, pl in sorted(geneset_plan.items()):
                sample = ", ".join(f"{k}->{v}" for k, v in list(pl["regulons"].items())[:2])
                _log(f"    {key:44} {len(pl['regulons'])} regulon names, "
                     f"{len(pl['targets'])} distinct target symbols "
                     f"(over {pl['n_targets']} entries)  ({sample})")
                if pl["unresolved"]:
                    _log(f"        {pl['unresolved']} names absent from var — left as-is")

        obsm_plan, obsm_skipped = plan_obsm_restore(args.adata, restored)
        if obsm_plan or obsm_skipped:
            _log("regulon-score column names in obsm:")
            for key, mapping in sorted(obsm_plan.items()):
                sample = ", ".join(f"{k}->{v}" for k, v in list(mapping.items())[:2])
                _log(f"    {key:44} {len(mapping):3} columns  ({sample})")
            for key, n, why in obsm_skipped:
                _log(f"    {key:44} {n:3} columns  LEFT AS-IS ({why})")
    else:
        geneset_plan = {}
        _log("--skip-uns: leaving uns tables and geneset dicts untouched")

    if args.dry_run:
        _log("dry run, nothing written")
        return 0
    if changed == 0 and not uns_plan and not geneset_plan and not obsm_plan:
        _log("nothing to change")
        return 0

    import scanpy as sc

    _log(f"reading {args.adata}")
    adata = sc.read_h5ad(args.adata)
    if list(adata.var_names) != target:
        raise SystemExit("var index changed between the check and the load — aborting")

    adata.var_names = restored

    for key, mapping in sorted(uns_plan.items()):
        table = adata.uns.get(key)
        if table is None or not hasattr(table, "rename"):
            _log(f"    uns['{key}'] is not a DataFrame in memory — skipped")
            continue
        before = list(table.index)
        table = table.rename(index=mapping)
        if len(set(table.index)) != len(before):
            # Renaming must never merge rows. Guaranteed by the collision check, but the
            # cost of being wrong here is silently dropped scores, so verify anyway.
            raise SystemExit(f"uns['{key}']: rename collapsed the index — aborting")
        adata.uns[key] = table
        _log(f"    uns['{key}']: renamed {len(mapping)} of {len(before)} rows")

    for key, pl in sorted(geneset_plan.items()):
        gs = adata.uns.get(key)
        if not isinstance(gs, dict):
            _log(f"    uns['{key}'] is not a dict in memory — skipped")
            continue
        rebuilt = {}
        for reg, tgts in gs.items():
            new_reg = pl["regulons"].get(reg, reg)
            if new_reg in rebuilt:
                raise SystemExit(f"uns['{key}']: regulon rename collided on {new_reg}")
            rebuilt[new_reg] = [pl["targets"].get(t, t) for t in list(tgts)]
        adata.uns[key] = rebuilt
        _log(f"    uns['{key}']: renamed {len(pl['regulons'])} regulons and "
             f"{len(pl['targets'])} distinct target symbols")

    for key, mapping in sorted(obsm_plan.items()):
        table = adata.obsm.get(key)
        if table is None or not hasattr(table, "rename"):
            _log(f"    obsm['{key}'] is not a DataFrame in memory — skipped")
            continue
        n_before = table.shape[1]
        table = table.rename(columns=mapping)
        if len(set(table.columns)) != n_before:
            raise SystemExit(f"obsm['{key}']: rename collapsed columns — aborting")
        adata.obsm[key] = table
        _log(f"    obsm['{key}']: renamed {len(mapping)} of {n_before} columns")

    # Keep a record on the object itself, so a later reader can tell this was repaired
    # rather than never broken.
    adata.uns["var_symbols_restored_from"] = os.path.basename(args.source)

    out_path = args.out or args.adata
    if os.path.abspath(out_path) == os.path.abspath(args.adata) and args.backup:
        backup = args.adata.replace(".h5ad", ".pre-symbol-restore.h5ad")
        if os.path.exists(backup):
            _log(f"backup already exists, leaving it: {backup}")
        else:
            _log(f"backing up to {backup}")
            shutil.copy2(args.adata, backup)

    # Same coercion the pipeline applies before writing.
    for key, value in adata.obsm.items():
        if hasattr(value, "columns"):
            value.columns = value.columns.astype(str)

    _log(f"writing {out_path}")
    adata.write(out_path)
    _log("done — regenerate the GeoJSON afterwards so the frontend sees the restored symbols")
    return 0


if __name__ == "__main__":
    sys.exit(main())
