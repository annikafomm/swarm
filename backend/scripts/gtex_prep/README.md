# Building GENIE3 / SPONGE networks from GTEx heart data

SWARM's upload form consumes **already-inferred regulatory networks**, not expression data.
It never runs GENIE3 or SPONGE itself — it only scores spatial data *against* a network you
supply. These scripts produce those networks for human heart tissue.

Without one of these networks, the entire **Regulatory Scores** tab (AUCell, GSVA, ssGSEA,
VIPER) and the Global Regulatory Scores table stay unavailable. That is the only score family
that cannot be backfilled from the spatial data alone.

## Pipeline

```
fetch_gtex_heart.py   download GTEx v10 + do the two tricky joins
        |
        +--> expr_symbols_<tissue>.csv ............ run_genie3.R --> genie3 network csv
        |
        +--> expr_ensembl_<tissue>_paired.csv .....\
        +--> expr_mirna_<tissue>_paired_mimat.csv ..} run_sponge.R --> 2 sponge csvs
        +--> sample_pairing_<tissue>.csv ........../

make_tf_list.py       CollecTRI TF universe, used as GENIE3's regulator set
```

## Why GTEx v10

It is the only openly available human-heart source that provides bulk RNA-seq **and**
sample-matched small RNA-seq. SPONGE's `mscor` is a correlation computed *across samples*, so
the mRNA and miRNA matrices must come from the same specimens — two separate cohorts of the
same tissue would not be statistically valid.

| Tissue | bulk RNA-seq | small RNA-seq | paired specimens |
|---|---|---|---|
| Heart – Left Ventricle | 452 | 399 | **361** |
| Heart – Atrial Appendage | 461 | 371 | **358** |

Fully open — no dbGaP application (only individual-level genotypes are controlled). Donors are
non-diseased postmortem, so it works as a healthy reference.

Small RNA exists **only in v10** — v11 has bulk RNA-seq only.

If your spatial data covers atrial or conduction-system regions, build a second network from
`--tissue atrial_appendage` rather than scoring everything against ventricle. There is no open
bulk cohort of human AVN/SAN conduction tissue at usable n, and GTEx has no conduction-system
tissue at all, so atrial appendage is the closest available proxy.

## The four joins that silently produce nothing if done wrong

Each of these was hit and fixed while building this; the scripts now fail loudly instead.

### 1. Identifier type differs per tool — symbols vs Ensembl

`backend/calc_R_scores/calc_scores.R` matches the two network types against different things:

| Network | Matched against | So ids must be |
|---|---|---|
| GENIE3 | `rownames(expr) <- var_df[,V1]` (AnnData var index) | **gene symbols** |
| SPONGE | `var_df[, get(ensembl_col)]` | **Ensembl gene ids** |

This is why `uns["genie_genesets"]` is keyed by symbols (`CEBPA`) while
`uns["sponge_genesets"]` is keyed by Ensembl (`ENSG00000034239`) — an asymmetry that is
otherwise baffling. Get it backwards and the pipeline logs
*"With the chosen parameters no regulons could be created."* rather than an error.

`fetch_gtex_heart.py` writes both flavours from the GCT's `Description` and `Name` columns.

### 2. GTEx miRNA rows are RNAcentral URS accessions

The miRNA matrix is keyed by `URS000012A1DD`, not by anything a miRNA-target database knows.
The annotation file (`smallRNA.filtered_annotated_031725.txt`) maps them.

### 3. TargetScan is keyed by MIMAT, not miRBase names

The non-obvious one. SPONGE's bundled `targetscan_ensg` has **348 columns, none matching
`^hsa-`** — they are MIMAT accessions (`MIMAT0000062`). Mapping URS to the human-readable
`hsa-miR-200b-5p` gives a **0-miRNA overlap** and SPONGE finds nothing.

Hence `--mirna-id mimat` is the default. Use `--mirna-id mirbase` only for eyeballing data.

### 4. mRNA and miRNA sample ids differ — pair on the specimen

They are separate *aliquots* of one tissue specimen:

```
GTEX-111YS-0426-SM-5987O   <- bulk RNA-seq aliquot
GTEX-111YS-0426-SM-G5KEX   <- small RNA-seq aliquot
GTEX-111YS-0426            <- the specimen, i.e. the join key
```

Intersecting the full sample ids yields **zero** matches. `fetch_gtex_heart.py` writes
`sample_pairing_<tissue>.csv`, and `run_sponge.R` uses it to align and reorder both matrices —
required, because `mscor` compares columns positionally.

## Usage

```bash
cd backend/scripts/gtex_prep

# 0. see what will be downloaded (~130 MB, cached afterwards)
python3 fetch_gtex_heart.py --tissue left_ventricle --for both --dry-run

# 1. download + prepare all matrices
python3 fetch_gtex_heart.py --tissue left_ventricle --for both --out-dir gtex_lv

# 2. GENIE3 regulator set from CollecTRI (same resource LIANA+ falls back to)
python3 make_tf_list.py --out gtex_lv/tfs.txt \
    --restrict-to gtex_lv/expr_symbols_left_ventricle.csv

# 3a. GENIE3 -> regulatoryGene,targetGene,weight
Rscript run_genie3.R --expr gtex_lv/expr_symbols_left_ventricle.csv \
    --tf-list gtex_lv/tfs.txt --max-targets 5000 --n-cores 32 \
    --out gtex_lv/genie3_heart_lv_network.csv

# 3b. SPONGE -> interactionNetwork.csv + networkAnalysis.csv
Rscript run_sponge.R --gene-expr gtex_lv/expr_ensembl_left_ventricle_paired.csv \
    --mir-expr gtex_lv/expr_mirna_left_ventricle_paired_mimat.csv \
    --pairing gtex_lv/sample_pairing_left_ventricle.csv \
    --out-dir gtex_lv --prefix sponge_heart_lv --n-cores 24
```

Both R scripts support `--dry-run`-style previews (`run_genie3.R --dry-run`) and print the
planned work before starting, which is worth checking since both are long-running.

## Cost control

**GENIE3** fits one random forest per target gene with the regulators as predictors, so it is
roughly `O(n_targets * n_regulators * n_samples)`. Handing it all ~58k GTEx genes as both
targets and regulators is not viable. Defaults:

- expression filter: >= 1 TPM in >= 20% of samples (57,853 -> ~14,400 genes on LV)
- regulators: the CollecTRI TF list, intersected with what survives the filter (~794)
- targets: 5,000 most variable, plus the regulators so TF-TF edges stay discoverable
- output truncated to the top 100,000 edges by weight, matching the pipeline's own `--top_n`
  default, which would truncate it anyway

**SPONGE** cost grows steeply with gene count too; `--max-genes 3000` keeps the most variable
genes. Both scripts register a `doParallel` backend — SPONGE parallelises via `foreach` and
runs single-threaded if no backend is registered, so `--n-cores` matters.

## Where the upload form wants each file

| File | Upload field |
|---|---|
| `genie3_*_network.csv` | Network Scores -> Upload Genie3 network |
| `sponge_*_interactionNetwork.csv` | Network Scores -> SPONGE "networkinteractions" |
| `sponge_*_networkAnalysis.csv` | Network Scores -> SPONGE "networkanalysis" |

GENIE3 unlocks VIPER plus AUCell/GSVA/ssGSEA; the SPONGE pair unlocks AUCell/GSVA/ssGSEA
(there is no VIPER for SPONGE). Both SPONGE files are required together — the form checks for
both before enabling submit.

## Requirements

- `SPONGE` (bundled with the image, 1.28.0) — ships TargetScan and miRcode, so no extra
  downloads
- `GENIE3` — **not** in the image. Install into a user library, since the conda R library is
  read-only:
  ```r
  dir.create("~/R/library", recursive = TRUE)
  BiocManager::install("GENIE3", lib = "~/R/library")
  ```
  `run_genie3.R` prepends `~/R/library` to `.libPaths()` itself.
- `decoupler` (Python, already installed) for `make_tf_list.py`; needs network access on first
  call to fetch CollecTRI.

## Fallback / cross-check cohorts

- **MAGNet, GSE141910** — 366 human LV (166 non-failing), open on GEO with ready NCBI
  counts/TPM TSVs. Good as an independent GENIE3 replicate to check network stability. Note
  the genotypes are controlled (dbGaP phs001539) but the expression data is not. No miRNA.
- **GSE46224** (n=40) and **GSE135055** (n=30) — the only other open, genuinely
  same-subject paired mRNA+miRNA human heart sets. At that n, `mscor` is underpowered; treat
  them as sanity checks, not replacements.
- The Heart Cell Atlas project itself (Kanemaru 2023, Litviňuková 2020) released **no bulk and
  no miRNA data** — only Visium/Multiome/single-cell. There is no same-project option.
