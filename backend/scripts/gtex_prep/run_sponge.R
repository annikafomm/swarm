#!/usr/bin/env Rscript
#
# Infer a SPONGE ceRNA network from paired gene + miRNA expression, and write the exact two
# CSVs SWARM's upload form expects.
#
# Why this exists
# ---------------
# SWARM does not run SPONGE — it consumes a pre-built network. `backend/calc_R_scores/calc_scores.R`
# validates the two uploads against a fixed column contract and will refuse anything else:
#
#   --sponge_network   geneA, geneB, df, cor, pcor, mscor, p.val, p.adj
#   --sponge_analysis  gene, degree, eigenvector, betweenness, page_rank
#
# (see the `interaction_cols` / `centrality_cols` check in calc_scores.R). This script emits
# precisely those, in that order.
#
# Identifier type matters and is easy to get wrong: calc_scores.R matches a SPONGE network
# against `var_df[, get(ensembl_col)]`, i.e. **Ensembl gene IDs** — not symbols. So geneA/geneB
# must be Ensembl. (GENIE3 is the opposite: it matches the AnnData var index, i.e. symbols.)
# `fetch_gtex_heart.py` already writes the gene matrix with Ensembl row ids for this reason.
#
# miRNA-target interactions come from SPONGE's own bundled TargetScan/miRcode tables, so no
# extra download is needed. The bundled tables come in both `_ensg` and `_symbol` flavours;
# we use `_ensg` to stay consistent with the Ensembl-keyed expression matrix.
#
# Usage
#   Rscript run_sponge.R --gene-expr expr_ensembl_left_ventricle_paired.csv \
#                        --mir-expr  expr_mirna_left_ventricle_paired.csv \
#                        --pairing   sample_pairing_left_ventricle.csv \
#                        --out-dir   . [--n-cores 8] [--max-genes 3000]

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

opts <- parse_args(OptionParser(option_list = list(
  make_option("--gene-expr", type = "character", dest = "gene_expr",
              help = "Gene expression csv, rows = Ensembl ids, cols = samples."),
  make_option("--mir-expr", type = "character", dest = "mir_expr",
              help = "miRNA expression csv, rows = miRBase names, cols = samples."),
  make_option("--pairing", type = "character", default = NULL,
              help = "specimen,bulk_sampid,mirna_sampid map used to align the two matrices."),
  make_option("--out-dir", type = "character", dest = "out_dir", default = ".",
              help = "Where to write the two csvs."),
  make_option("--prefix", type = "character", default = "sponge",
              help = "Output filename prefix [default %default]."),
  make_option("--n-cores", type = "integer", dest = "n_cores", default = 4L,
              help = "Cores for SPONGE's parallel steps [default %default]."),
  make_option("--max-genes", type = "integer", dest = "max_genes", default = 3000L,
              help = paste("Keep only the N most variable genes. SPONGE cost grows steeply",
                           "with gene count; 0 disables the cap [default %default].")),
  make_option("--mscor-cutoff", type = "double", dest = "mscor_cutoff", default = 0.1,
              help = "Report interactions with mscor above this [default %default]."),
  make_option("--null-perms", type = "integer", dest = "null_perms", default = 1e5,
              help = "Null-model samples for p-values [default %default].")
)))

for (req in c("gene_expr", "mir_expr")) {
  if (is.null(opts[[req]])) stop(sprintf("--%s is required", gsub("_", "-", req)))
}

log_msg <- function(...) cat(sprintf("[sponge] %s\n", paste0(...)), sep = "")
t_start <- Sys.time()
took <- function(t0) sprintf("%.1fs", as.numeric(difftime(Sys.time(), t0, units = "secs")))

suppressPackageStartupMessages({
  library(SPONGE)
  library(igraph)
})

# SPONGE parallelises through foreach, which needs a backend registered by the caller —
# without this every step runs single-threaded regardless of --n-cores.
if (opts$n_cores > 1L) {
  ok <- requireNamespace("doParallel", quietly = TRUE)
  if (ok) {
    cl <- parallel::makeCluster(opts$n_cores)
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    cat(sprintf("[sponge] registered %d-core parallel backend\n", opts$n_cores))
  } else {
    cat("[sponge] doParallel not available; running single-threaded\n")
  }
}

# ---- load and align -----------------------------------------------------------------

log_msg("reading expression matrices ...")
gene_dt <- fread(opts$gene_expr)
mir_dt  <- fread(opts$mir_expr)

gene_ids <- gene_dt[[1]]
mir_ids  <- mir_dt[[1]]
gene_mat <- as.matrix(gene_dt[, -1, with = FALSE])
mir_mat  <- as.matrix(mir_dt[, -1, with = FALSE])
rownames(gene_mat) <- gene_ids
rownames(mir_mat)  <- mir_ids
log_msg(sprintf("  genes: %d x %d | miRNAs: %d x %d",
                nrow(gene_mat), ncol(gene_mat), nrow(mir_mat), ncol(mir_mat)))

# SPONGE expects samples in rows, and — critically — the same samples in the same order in
# both matrices, since mscor is a correlation across samples. The two GTEx matrices are keyed
# by different aliquot ids, so align via the specimen pairing map rather than by name.
if (!is.null(opts$pairing)) {
  pair <- fread(opts$pairing)
  keep <- pair[bulk_sampid %in% colnames(gene_mat) & mirna_sampid %in% colnames(mir_mat)]
  log_msg(sprintf("  pairing map: %d of %d specimens present in both matrices",
                  nrow(keep), nrow(pair)))
  if (nrow(keep) < 20) {
    stop("fewer than 20 paired specimens — mscor would be meaningless at this n")
  }
  gene_mat <- gene_mat[, keep$bulk_sampid, drop = FALSE]
  mir_mat  <- mir_mat[,  keep$mirna_sampid, drop = FALSE]
  # Rename both to the neutral specimen id so the transposed matrices line up exactly.
  colnames(gene_mat) <- keep$specimen
  colnames(mir_mat)  <- keep$specimen
} else {
  common <- intersect(colnames(gene_mat), colnames(mir_mat))
  if (length(common) < 20) {
    stop(paste("no --pairing given and only", length(common),
               "column names are shared; GTEx aliquot ids differ between assays,",
               "so pass --pairing"))
  }
  gene_mat <- gene_mat[, common, drop = FALSE]
  mir_mat  <- mir_mat[,  common, drop = FALSE]
}
stopifnot(identical(colnames(gene_mat), colnames(mir_mat)))
log_msg(sprintf("  aligned on %d specimens", ncol(gene_mat)))

# ---- filter ------------------------------------------------------------------------

# log-transform TPM: SPONGE's partial correlations assume roughly linear relationships, and
# raw TPM is heavily skewed.
gene_mat <- log2(gene_mat + 1)
mir_mat  <- log2(mir_mat + 1)

drop_flat <- function(m, label) {
  keep <- apply(m, 1, function(r) sd(r, na.rm = TRUE) > 0)
  if (any(!keep)) log_msg(sprintf("  dropped %d zero-variance %s", sum(!keep), label))
  m[keep, , drop = FALSE]
}
gene_mat <- drop_flat(gene_mat, "genes")
mir_mat  <- drop_flat(mir_mat,  "miRNAs")

if (opts$max_genes > 0 && nrow(gene_mat) > opts$max_genes) {
  vars <- apply(gene_mat, 1, var)
  gene_mat <- gene_mat[order(vars, decreasing = TRUE)[seq_len(opts$max_genes)], , drop = FALSE]
  log_msg(sprintf("  kept %d most variable genes (--max-genes)", nrow(gene_mat)))
}

# SPONGE wants samples x features.
gene_expr <- t(gene_mat)
mir_expr  <- t(mir_mat)

# ---- miRNA-target candidates -------------------------------------------------------

# Bundled TargetScan, Ensembl-keyed on rows to match our gene ids. SPONGE exports this
# lazily, so reference it directly rather than via data() — data() emits a spurious
# "data set not found" warning even though the object resolves fine.
mir_predicted <- SPONGE::targetscan_ensg
log_msg(sprintf("bundled TargetScan (ensg): %d genes x %d miRNAs",
                nrow(mir_predicted), ncol(mir_predicted)))

overlap_g <- length(intersect(colnames(gene_expr), rownames(mir_predicted)))
overlap_m <- length(intersect(colnames(mir_expr),  colnames(mir_predicted)))
log_msg(sprintf("  overlap with our data: %d genes, %d miRNAs", overlap_g, overlap_m))
if (overlap_g == 0 || overlap_m == 0) {
  # The miRNA half is the usual culprit: TargetScan's columns are MIMAT accessions
  # (MIMAT0000062), NOT miRBase names (hsa-miR-200b-5p). Verified: 0 of its 348 columns
  # match ^hsa-. fetch_gtex_heart.py must therefore be run with --mirna-id mimat.
  stop(sprintf(paste(
    "no overlap with TargetScan (genes: %d, miRNAs: %d).",
    "Genes must be Ensembl ids (ENSG...) — ours look like '%s'.",
    "miRNAs must be MIMAT accessions (MIMAT...) — ours look like '%s'.",
    "If the miRNA ids are hsa-miR-* names, re-run fetch_gtex_heart.py with --mirna-id mimat."),
    overlap_g, overlap_m, colnames(gene_expr)[1], colnames(mir_expr)[1]))
}

t0 <- Sys.time()
log_msg("filtering gene-miRNA candidate interactions (elastic net) ...")
genes_mir_candidates <- sponge_gene_miRNA_interaction_filter(
  gene_expr = gene_expr,
  mir_expr = mir_expr,
  mir_predicted_targets = mir_predicted
)
log_msg(sprintf("  done in %s", took(t0)))

# ---- ceRNA interactions ------------------------------------------------------------

t0 <- Sys.time()
log_msg("computing ceRNA interactions (this is the expensive step) ...")
ceRNA <- sponge(
  gene_expr = gene_expr,
  mir_expr = mir_expr,
  mir_interactions = genes_mir_candidates
)
log_msg(sprintf("  %d candidate interactions in %s", nrow(ceRNA), took(t0)))
if (nrow(ceRNA) == 0) stop("SPONGE returned no interactions; nothing to write")

t0 <- Sys.time()
log_msg("building null model and computing p-values ...")
null_model <- sponge_build_null_model(
  number_of_datasets = opts$null_perms,
  number_of_samples = nrow(gene_expr)
)
ceRNA_p <- sponge_compute_p_values(sponge_result = ceRNA, null_model = null_model)
log_msg(sprintf("  done in %s", took(t0)))

res <- as.data.table(ceRNA_p)

# calc_scores.R requires exactly these names. SPONGE's own output already uses geneA/geneB/df/
# cor/pcor/mscor; the p-value columns vary by version, so normalise them defensively.
if (!"p.val" %in% names(res)) {
  cand <- intersect(c("p.value", "pval", "p_val"), names(res))
  if (length(cand)) setnames(res, cand[1], "p.val")
}
if (!"p.adj" %in% names(res)) {
  cand <- intersect(c("p.adj.value", "padj", "p_adj", "p.adjusted"), names(res))
  if (length(cand)) setnames(res, cand[1], "p.adj")
  else res[, p.adj := p.adjust(p.val, method = "BH")]
}

required <- c("geneA", "geneB", "df", "cor", "pcor", "mscor", "p.val", "p.adj")
missing <- setdiff(required, names(res))
if (length(missing)) {
  stop(sprintf("SPONGE output lacks required column(s): %s (present: %s)",
               paste(missing, collapse = ", "), paste(names(res), collapse = ", ")))
}

interactions <- res[, ..required]
n_before <- nrow(interactions)
interactions <- interactions[!is.na(mscor) & mscor > opts$mscor_cutoff]
log_msg(sprintf("kept %d of %d interactions with mscor > %s",
                nrow(interactions), n_before, opts$mscor_cutoff))
if (nrow(interactions) == 0) stop("no interactions passed the mscor cutoff; lower --mscor-cutoff")

# ---- centralities ------------------------------------------------------------------

t0 <- Sys.time()
log_msg("computing node centralities ...")
# sponge_node_centralities returns degree/eigenvector/betweenness/page_rank per gene, which is
# exactly the centrality contract. Feed it the filtered set so the two files describe the
# same network.
centralities <- as.data.table(sponge_node_centralities(interactions))
log_msg(sprintf("  %d nodes in %s", nrow(centralities), took(t0)))

if (!"gene" %in% names(centralities)) {
  cand <- intersect(c("geneA", "node", "name", "V1"), names(centralities))
  if (length(cand)) setnames(centralities, cand[1], "gene")
  else centralities[, gene := rownames(centralities)]
}
cent_required <- c("gene", "degree", "eigenvector", "betweenness", "page_rank")
cent_missing <- setdiff(cent_required, names(centralities))
if (length(cent_missing)) {
  stop(sprintf("centralities lack required column(s): %s (present: %s)",
               paste(cent_missing, collapse = ", "), paste(names(centralities), collapse = ", ")))
}
centralities <- centralities[, ..cent_required]

# ---- write ------------------------------------------------------------------------

dir.create(opts$out_dir, showWarnings = FALSE, recursive = TRUE)
f_int  <- file.path(opts$out_dir, sprintf("%s_interactionNetwork.csv", opts$prefix))
f_cent <- file.path(opts$out_dir, sprintf("%s_networkAnalysis.csv", opts$prefix))
fwrite(interactions, f_int)
fwrite(centralities, f_cent)

log_msg(sprintf("wrote %s (%d rows)", f_int, nrow(interactions)))
log_msg(sprintf("wrote %s (%d rows)", f_cent, nrow(centralities)))
log_msg("upload these as --sponge_network and --sponge_analysis respectively")
log_msg(sprintf("total runtime %s", took(t_start)))
