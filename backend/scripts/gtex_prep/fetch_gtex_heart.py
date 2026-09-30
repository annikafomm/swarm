#!/usr/bin/env python3
"""
Download GTEx v10 heart expression data and prepare matrices for GENIE3 / SPONGE.

Why this exists
---------------
SWARM's upload form takes *already-inferred networks*, not expression data: a GENIE3
`regulatoryGene,targetGene,weight` csv, and/or a SPONGE pair of
`geneA,geneB,df,cor,pcor,mscor,p.val,p.adj` + `gene,degree,eigenvector,betweenness,page_rank`.
Producing those for a new tissue means (1) getting a bulk expression matrix, (2) for SPONGE
also a *sample-matched* miRNA matrix, then (3) running the inference. This script covers
steps 1-2; `run_genie3.R` and `run_sponge.R` next to it cover step 3.

GTEx v10 is the only open human-heart source that serves both halves from the same tissue
specimens (verified sample counts: LV 452 bulk / 399 miRNA / **361 paired**; atrial appendage
461 / 371 / **358**). SPONGE's mscor is computed across samples, so the pairing has to be
genuine — two separate cohorts would not be valid.

Two joins this script exists to get right, both of which fail silently if done naively:

1. **Identifier type differs per tool.** The R pipeline matches a GENIE3 network against
   `rownames(expr) <- var_df[,V1]` (the AnnData var index, i.e. **gene symbols**), but matches
   a SPONGE network against `var_df[, get(ensembl_col)]` (**Ensembl gene IDs**) — see
   `backend/calc_R_scores/calc_scores.R`. GTEx's GCT carries both (`Name` = versioned Ensembl,
   `Description` = symbol), so this writes symbol-indexed output for GENIE3 and
   Ensembl-indexed (version stripped) output for SPONGE.
2. **GTEx miRNA rows are RNAcentral URS accessions, not miRBase names** (`URS000012A1DD`),
   and **miRNA/mRNA sample IDs differ** because they are separate aliquots of one specimen
   (`GTEX-1117F-0226-SM-DEURW` vs `...-SM-5GZWO`). Pairing is on the tissue-sample prefix
   (`GTEX-1117F-0226`); the URS -> `hsa-miR-*` map comes from the annotation file.

Downloads are cached in --cache-dir and skipped when the expected size already matches, so
re-runs are cheap.

Examples
--------
    # what would be downloaded, and how big
    python3 fetch_gtex_heart.py --tissue left_ventricle --dry-run

    # GENIE3 input only (no miRNA download -> much smaller)
    python3 fetch_gtex_heart.py --tissue left_ventricle --for genie3 --out-dir gtex_lv

    # both, including the paired miRNA matrix for SPONGE
    python3 fetch_gtex_heart.py --tissue left_ventricle --for both --out-dir gtex_lv
"""

import argparse
import gzip
import os
import sys
import urllib.request

GTEX_BASE = "https://storage.googleapis.com/adult-gtex"

TISSUES = {
    "left_ventricle": "heart_left_ventricle",
    "atrial_appendage": "heart_atrial_appendage",
}

# Sizes are the values observed on 2026-08-21; used only to decide whether a cached file
# looks complete, never to reject an otherwise-valid download.
FILES = {
    "bulk_tpm": {
        "url": f"{GTEX_BASE}/bulk-gex/v10/rna-seq/tpms-by-tissue/gene_tpm_v10_{{tissue}}.gct.gz",
        "approx_mb": {"heart_left_ventricle": 48.7, "heart_atrial_appendage": 52.9},
        "needed_for": ("genie3", "sponge"),
    },
    "mirna_tpm": {
        # One file for every tissue (16,761 samples); we subset the heart columns locally.
        "url": f"{GTEX_BASE}/bulk-gex/v10/small-RNA/miRNA_TPM_matrix_PORTAL_2025_03_17.txt.gz",
        "approx_mb": 61.6,
        "needed_for": ("sponge",),
    },
    "smallrna_annot": {
        # Maps RNAcentral URS ids -> miRBase names / MIMAT accessions.
        "url": f"{GTEX_BASE}/annotations/v10/small-RNA/smallRNA.filtered_annotated_031725.txt",
        "approx_mb": 7.9,
        "needed_for": ("sponge",),
    },
    "sample_attrs": {
        # SAMPID -> SMTSD (tissue) and SMAFRZE (RNASEQ / SMLRNA), needed to pick heart columns.
        "url": f"{GTEX_BASE}/annotations/v10/metadata-files/GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt",
        "approx_mb": 12.0,
        "needed_for": ("genie3", "sponge"),
    },
}


def _log(msg: str) -> None:
    print(f"[gtex-prep] {msg}", flush=True)


def _human(n_bytes: int) -> str:
    return f"{n_bytes / 1048576:.1f} MB"


def download(url: str, dest: str, expect_mb: float | None = None) -> str:
    """Fetch `url` to `dest`, skipping when a plausibly-complete copy is already there."""
    if os.path.exists(dest):
        have = os.path.getsize(dest)
        # Tolerate a wide band: the point is to catch truncated/aborted downloads, not to
        # pin an exact byte count that will drift when GTEx re-releases a file.
        if expect_mb is None or have > expect_mb * 1048576 * 0.5:
            _log(f"cached {os.path.basename(dest)} ({_human(have)})")
            return dest
        _log(f"cached copy looks truncated ({_human(have)}), re-downloading")

    tmp = dest + ".part"
    _log(f"downloading {os.path.basename(dest)} ...")

    def _progress(block_num, block_size, total_size):
        if total_size <= 0:
            return
        pct = min(100, block_num * block_size * 100 // total_size)
        if pct % 20 == 0 and block_num * block_size % (total_size // 5 or 1) < block_size:
            print(f"    {pct}%", end="\r", flush=True)

    urllib.request.urlretrieve(url, tmp, reporthook=_progress)
    os.replace(tmp, dest)
    _log(f"got {os.path.basename(dest)} ({_human(os.path.getsize(dest))})")
    return dest


def read_sample_attrs(path: str, tissue_smtsd: str):
    """Return (bulk_sampids, mirna_sampids) for one tissue.

    SMAFRZE marks which analysis freeze a sample belongs to: RNASEQ for bulk, SMLRNA for
    small RNA. Filtering on it is what makes the counts line up with the published numbers.
    """
    bulk, mirna = [], []
    with open(path, encoding="utf-8", errors="replace") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        i_id, i_tissue, i_frz = (
            header.index("SAMPID"),
            header.index("SMTSD"),
            header.index("SMAFRZE"),
        )
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= max(i_id, i_tissue, i_frz):
                continue
            if parts[i_tissue] != tissue_smtsd:
                continue
            if parts[i_frz] == "RNASEQ":
                bulk.append(parts[i_id])
            elif parts[i_frz] == "SMLRNA":
                mirna.append(parts[i_id])
    return bulk, mirna


def specimen_of(sampid: str) -> str:
    """GTEX-1117F-0226-SM-DEURW -> GTEX-1117F-0226.

    The trailing `-SM-<aliquot>` identifies the aliquot, not the tissue specimen. mRNA and
    miRNA were prepared from different aliquots of the same specimen, so the prefix is the
    only correct join key between the two matrices.
    """
    parts = sampid.split("-")
    return "-".join(parts[:3]) if len(parts) >= 3 else sampid


def read_urs_map(path: str, id_kind: str) -> dict:
    """Map RNAcentral URS id -> MIMAT accession or miRBase name, for miRNA rows only.

    Which one you want depends entirely on what the downstream network tool keys on, and
    getting it wrong produces a silent zero-overlap join:

      * **mimat** (default) — SPONGE's bundled TargetScan/miRcode tables are keyed by MIMAT
        accession (`MIMAT0000062`), *not* by miRBase name. Verified against
        `colnames(SPONGE::targetscan_ensg)`: 348 columns, 0 of them matching `^hsa-`.
      * **mirbase** — human-readable `hsa-miR-200b-5p`. Use for inspection, or for a tool
        that keys on names.
    """
    col = {"mimat": "MIMAT", "mirbase": "mirbase_name"}[id_kind]
    mapping = {}
    with open(path, encoding="utf-8", errors="replace") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        i_id = header.index("id") if "id" in header else 0
        if col not in header:
            raise SystemExit(
                f"no '{col}' column in {path}; columns were: {header[:12]}"
            )
        i_val = header.index(col)
        i_type = header.index("type") if "type" in header else None
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) <= max(i_id, i_val):
                continue
            if i_type is not None and len(p) > i_type and p[i_type] and p[i_type] != "miRNA":
                continue
            val = p[i_val].strip()
            if val and val.lower() not in ("na", "nan", "none", "-"):
                mapping[p[i_id].strip()] = val
    return mapping


def _open_text(path: str):
    return gzip.open(path, "rt", encoding="utf-8", errors="replace") if path.endswith(".gz") \
        else open(path, encoding="utf-8", errors="replace")


def subset_gct(gct_path: str, keep_samples: set, out_path: str, id_mode: str) -> tuple:
    """Stream a GTEx GCT to a plain csv of genes x kept samples.

    GCT layout: line 1 `#1.2`, line 2 `<n_rows>\t<n_cols>`, line 3 header
    (`Name`, `Description`, then sample ids), then one row per gene.

    id_mode:
      "symbol"  -> row id from `Description` (what a GENIE3 network must use)
      "ensembl" -> row id from `Name` with the version suffix stripped (what SPONGE must use)

    Streamed line-by-line rather than loaded into a DataFrame: these files are ~60k rows x
    ~460 columns and there is no need to hold them in memory to do a column subset.
    """
    kept_idx, out_names, n_rows = None, None, 0
    seen = set()
    with _open_text(gct_path) as fh, open(out_path, "w", encoding="utf-8") as out:
        fh.readline()  # "#1.2"
        fh.readline()  # dims
        header = fh.readline().rstrip("\n").split("\t")
        if header[:2] != ["Name", "Description"]:
            raise SystemExit(f"unexpected GCT header in {gct_path}: {header[:4]}")
        kept_idx = [i for i, s in enumerate(header) if i >= 2 and s in keep_samples]
        out_names = [header[i] for i in kept_idx]
        if not kept_idx:
            raise SystemExit(
                "no requested samples found in the GCT — check the tissue name matches "
                "SMTSD in the sample attributes file"
            )
        out.write("gene," + ",".join(out_names) + "\n")
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 3:
                continue
            gid = p[1].strip() if id_mode == "symbol" else p[0].split(".")[0]
            if not gid or gid in seen:
                # Duplicate symbols exist in GTEx (several Ensembl ids share a Description);
                # keeping the first occurrence keeps row ids unique, which both downstream
                # tools require.
                continue
            seen.add(gid)
            out.write(gid + "," + ",".join(p[i] for i in kept_idx) + "\n")
            n_rows += 1
    return n_rows, out_names


def subset_mirna(mirna_path: str, keep_samples: set, urs_map: dict, out_path: str) -> tuple:
    """Subset the all-tissue miRNA TPM matrix and rename URS ids to miRBase names."""
    n_rows, unmapped = 0, 0
    with _open_text(mirna_path) as fh, open(out_path, "w", encoding="utf-8") as out:
        header = fh.readline().rstrip("\n").split("\t")
        kept_idx = [i for i, s in enumerate(header) if i >= 1 and s in keep_samples]
        if not kept_idx:
            raise SystemExit("no requested miRNA samples found in the miRNA matrix")
        out.write("miRNA," + ",".join(header[i] for i in kept_idx) + "\n")
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 2:
                continue
            name = urs_map.get(p[0].strip())
            if not name:
                unmapped += 1
                continue
            out.write(name + "," + ",".join(p[i] for i in kept_idx) + "\n")
            n_rows += 1
    return n_rows, unmapped, [header[i] for i in kept_idx]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--tissue", default="left_ventricle", choices=sorted(TISSUES))
    ap.add_argument("--for", dest="target", default="both",
                    choices=["genie3", "sponge", "both"])
    ap.add_argument("--out-dir", default="gtex_heart")
    ap.add_argument("--cache-dir", default=None,
                    help="Where downloads live (default: <out-dir>/_downloads).")
    ap.add_argument("--mirna-id", default="mimat", choices=["mimat", "mirbase"],
                    help="Row id for the miRNA matrix. 'mimat' (default) is what SPONGE's "
                         "bundled TargetScan is keyed by; 'mirbase' is human-readable.")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    tissue_slug = TISSUES[args.tissue]
    smtsd = {"left_ventricle": "Heart - Left Ventricle",
             "atrial_appendage": "Heart - Atrial Appendage"}[args.tissue]
    wants = {"genie3": {"genie3"}, "sponge": {"sponge"}, "both": {"genie3", "sponge"}}[args.target]

    cache = args.cache_dir or os.path.join(args.out_dir, "_downloads")
    os.makedirs(cache, exist_ok=True)
    os.makedirs(args.out_dir, exist_ok=True)

    plan = []
    for key, spec in FILES.items():
        if not (set(spec["needed_for"]) & wants):
            continue
        url = spec["url"].format(tissue=tissue_slug)
        mb = spec["approx_mb"]
        if isinstance(mb, dict):
            mb = mb.get(tissue_slug)
        plan.append((key, url, os.path.join(cache, os.path.basename(url)), mb))

    _log(f"tissue: {smtsd}  |  target: {args.target}")
    total = sum(mb for *_, mb in plan if mb)
    for key, url, dest, mb in plan:
        state = "cached" if os.path.exists(dest) else "to download"
        _log(f"  {key:16} {mb or '?':>6} MB  {state}")
    _log(f"  total download if nothing cached: ~{total:.0f} MB")
    if args.dry_run:
        _log("dry run, nothing fetched")
        return 0

    paths = {key: download(url, dest, mb) for key, url, dest, mb in plan}

    _log("reading sample annotations ...")
    bulk_ids, mirna_ids = read_sample_attrs(paths["sample_attrs"], smtsd)
    _log(f"  {smtsd}: {len(bulk_ids)} bulk RNA-seq, {len(mirna_ids)} small RNA-seq samples")

    # Pair on tissue specimen, not aliquot.
    bulk_by_spec = {specimen_of(s): s for s in bulk_ids}
    mirna_by_spec = {specimen_of(s): s for s in mirna_ids}
    paired_spec = sorted(set(bulk_by_spec) & set(mirna_by_spec))
    _log(f"  paired specimens (same tissue sample, both assays): {len(paired_spec)}")

    written = {}

    if "genie3" in wants:
        # GENIE3 gets every bulk sample — more samples is strictly better for it, and it has
        # no need to match the miRNA side.
        out = os.path.join(args.out_dir, f"expr_symbols_{args.tissue}.csv")
        n, cols = subset_gct(paths["bulk_tpm"], set(bulk_ids), out, id_mode="symbol")
        _log(f"GENIE3 input: {out}  ({n} genes x {len(cols)} samples, gene symbols)")
        written["genie3_expr"] = out

    if "sponge" in wants:
        # SPONGE needs both matrices restricted to the *same* specimens, in the same order.
        keep_bulk = {bulk_by_spec[s] for s in paired_spec}
        keep_mirna = {mirna_by_spec[s] for s in paired_spec}

        out_g = os.path.join(args.out_dir, f"expr_ensembl_{args.tissue}_paired.csv")
        n_g, cols_g = subset_gct(paths["bulk_tpm"], keep_bulk, out_g, id_mode="ensembl")
        _log(f"SPONGE gene input: {out_g}  ({n_g} genes x {len(cols_g)} samples, Ensembl ids)")

        urs_map = read_urs_map(paths["smallrna_annot"], args.mirna_id)
        _log(f"  URS -> {args.mirna_id} map: {len(urs_map)} miRNA entries")
        out_m = os.path.join(
            args.out_dir, f"expr_mirna_{args.tissue}_paired_{args.mirna_id}.csv")
        n_m, unmapped, cols_m = subset_mirna(paths["mirna_tpm"], keep_mirna, urs_map, out_m)
        _log(f"SPONGE miRNA input: {out_m}  ({n_m} miRNAs x {len(cols_m)} samples, {args.mirna_id} ids)")
        if unmapped:
            _log(f"  dropped {unmapped} small-RNA rows with no miRBase name (piRNA/snRNA/etc.)")

        # The two matrices must be column-aligned by specimen for mscor to be meaningful.
        keymap = os.path.join(args.out_dir, f"sample_pairing_{args.tissue}.csv")
        with open(keymap, "w", encoding="utf-8") as fh:
            fh.write("specimen,bulk_sampid,mirna_sampid\n")
            for s in paired_spec:
                fh.write(f"{s},{bulk_by_spec[s]},{mirna_by_spec[s]}\n")
        _log(f"sample pairing map: {keymap}  ({len(paired_spec)} rows)")
        written.update({"sponge_gene": out_g, "sponge_mirna": out_m, "pairing": keymap})

    _log("done. next step:")
    if "genie3" in wants:
        _log(f"    Rscript run_genie3.R --expr {written['genie3_expr']} --out genie3_network.csv")
    if "sponge" in wants:
        _log(f"    Rscript run_sponge.R --gene-expr {written['sponge_gene']} "
             f"--mir-expr {written['sponge_mirna']} --pairing {written['pairing']} --out-dir .")
    return 0


if __name__ == "__main__":
    sys.exit(main())
