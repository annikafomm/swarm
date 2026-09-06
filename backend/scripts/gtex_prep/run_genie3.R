#!/usr/bin/env Rscript
#
# Infer a GENIE3 regulatory network from a bulk expression matrix and write the exact CSV
# SWARM's upload form expects.
#
# Why this exists
# ---------------
# SWARM consumes a pre-built GENIE3 network rather than inferring one. The contract is set by
# `backend/calc_R_scores/utils.R:create_Genie_modules`, which renames columns by position-free
# name lookup:
#
#     setnames(regulon_df, c("regulatoryGene", "targetGene", "weight"), c("tf", "target", "mi"))
#
# so the upload must have exactly those three column names. This script emits them.
#
# **Identifiers must be gene SYMBOLS.** calc_scores.R sets `rownames(expr) <- var_df[,V1]`
# (the AnnData var index) before matching the network, and that index holds symbols. This is
# the opposite of the SPONGE path, which matches on Ensembl ids — see run_sponge.R. Feeding a
# symbol network where Ensembl is expected (or vice versa) yields zero regulons and a
# "With the chosen parameters no regulons could be created." message rather than an error.
#
# Cost: GENIE3 fits one random forest per target gene with the regulators as predictors, so
# runtime is roughly O(n_targets * n_regulators * n_samples). On 452 GTEx samples, restricting
# regulators to a TF list (see make_tf_list.py) and targets to the most variable genes is what
# makes this finish in minutes rather than days. Both are capped by default here.
#
# Usage
#   Rscript run_genie3.R --expr expr_symbols_left_ventricle.csv --out genie3_network.csv \
#       [--tf-list tfs.txt] [--max-targets 5000] [--n-cores 8] [--top-n 100000]

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

opts <- parse_args(OptionParser(option_list = list(
  make_option("--expr", type = "character",
              help = "Expression csv: first column gene SYMBOL, remaining columns samples."),
  make_option("--out", type = "character", default = "genie3_network.csv",
              help = "Output network csv [default %default]."),
  make_option("--tf-list", type = "character", dest = "tf_list", default = NULL,
              help = paste("File with one TF symbol per line, used as GENIE3 regulators.",
                           "Strongly recommended; without it every gene is a candidate",
                           "regulator and the run becomes impractical.")),
  make_option("--max-targets", type = "integer", dest = "max_targets", default = 5000L,
              help = "Keep the N most variable genes as targets; 0 = all [default %default]."),
  make_option("--min-tpm", type = "double", dest = "min_tpm", default = 1.0,
              help = "Drop genes not reaching this TPM in at least --min-frac of samples."),
  make_option("--min-frac", type = "double", dest = "min_frac", default = 0.2,
              help = "Fraction of samples for the --min-tpm filter [default %default]."),
  make_option("--top-n", type = "integer", dest = "top_n", default = 100000L,
              help = paste("Write only the top N edges by weight. Matches the pipeline's own",
                           "--top_n default, which truncates the network anyway [default %default].")),
  make_option("--n-cores", type = "integer", dest = "n_cores", default = 4L,
              help = "Cores for GENIE3 [default %default]."),
  make_option("--n-trees", type = "integer", dest = "n_trees", default = 500L,
              help = paste("Trees per random forest. GENIE3's own default is 1000;",
                           "500 roughly halves runtime with little ranking change",
                           "[default %default].")),
  make_option("--seed", type = "integer", default = 42L, help = "RNG seed [default %default]."),
  make_option("--dry-run", action = "store_true", dest = "dry_run", default = FALSE,
              help = "Report the matrix shape and planned work, then exit.")
)))

if (is.null(opts$expr)) stop("--expr is required")

log_msg <- function(...) cat(sprintf("[genie3] %s\n", paste0(...)), sep = "")
t_start <- Sys.time()
took <- function(t0) sprintf("%.1fs", as.numeric(difftime(Sys.time(), t0, units = "secs")))

# GENIE3 was installed into a user library, which non-interactive Rscript does not always pick
# up depending on how R_LIBS_USER is expanded.
user_lib <- path.expand("~/R/library")
if (dir.exists(user_lib)) .libPaths(c(user_lib, .libPaths()))

log_msg("reading expression matrix ...")
t0 <- Sys.time()
dt <- fread(opts$expr)
gene_ids <- dt[[1]]
mat <- as.matrix(dt[, -1, with = FALSE])
rownames(mat) <- gene_ids
log_msg(sprintf("  %d genes x %d samples in %s", nrow(mat), ncol(mat), took(t0)))

# ---- filtering ---------------------------------------------------------------------

# Expression filter first: GENIE3 on near-zero rows produces noise-driven edges.
if (opts$min_tpm > 0) {
  keep <- rowSums(mat >= opts$min_tpm) >= (opts$min_frac * ncol(mat))
  log_msg(sprintf("  expression filter (>= %g TPM in >= %.0f%% of samples): %d of %d genes kept",
                  opts$min_tpm, opts$min_frac * 100, sum(keep), nrow(mat)))
  mat <- mat[keep, , drop = FALSE]
}

# log-transform: GENIE3 is tree-based so monotone transforms do not change splits, but it
# keeps variance ranking below from being dominated by a few very high-TPM genes.
mat <- log2(mat + 1)

keep_var <- apply(mat, 1, function(r) sd(r) > 0)
if (any(!keep_var)) {
  log_msg(sprintf("  dropped %d zero-variance genes", sum(!keep_var)))
  mat <- mat[keep_var, , drop = FALSE]
}

regulators <- NULL
if (!is.null(opts$tf_list)) {
  tfs <- readLines(opts$tf_list, warn = FALSE)
  tfs <- unique(trimws(tfs[nzchar(trimws(tfs))]))
  regulators <- intersect(tfs, rownames(mat))
  log_msg(sprintf("  regulators: %d of %d listed TFs present in the matrix",
                  length(regulators), length(tfs)))
  if (length(regulators) < 10) {
    stop(paste("fewer than 10 TFs matched the matrix rows. The matrix must be keyed by gene",
               "SYMBOL for this to work — check you passed the symbol matrix, not Ensembl."))
  }
} else {
  log_msg("  no --tf-list given: every gene is a candidate regulator (slow)")
}

# Cap targets by variance, but always keep the regulators themselves as candidate targets so
# TF-TF edges remain discoverable.
if (opts$max_targets > 0 && nrow(mat) > opts$max_targets) {
  vars <- apply(mat, 1, var)
  top <- names(sort(vars, decreasing = TRUE))[seq_len(opts$max_targets)]
  keep_rows <- union(top, if (is.null(regulators)) character(0) else regulators)
  mat <- mat[intersect(rownames(mat), keep_rows), , drop = FALSE]
  log_msg(sprintf("  targets capped to %d most variable (+ regulators) -> %d rows",
                  opts$max_targets, nrow(mat)))
}
if (!is.null(regulators)) regulators <- intersect(regulators, rownames(mat))

n_reg <- if (is.null(regulators)) nrow(mat) else length(regulators)
log_msg(sprintf("planned: %d targets x %d regulators x %d samples, %d core(s)",
                nrow(mat), n_reg, ncol(mat), opts$n_cores))
if (opts$dry_run) { log_msg("dry run, stopping"); quit(status = 0) }

# ---- inference --------------------------------------------------------------------

suppressPackageStartupMessages(library(GENIE3))

t0 <- Sys.time()
log_msg("running GENIE3 (this is the long step) ...")
# GENIE3 1.28.0 takes no `seed` argument (signature: exprMatrix, regulators, targets,
# treeMethod, K, nTrees, nCores, returnMatrix, verbose) — seed externally instead.
set.seed(opts$seed)
wm <- GENIE3(
  exprMatrix = mat,
  regulators = regulators,
  nTrees = opts$n_trees,
  nCores = opts$n_cores,
  verbose = TRUE
)
log_msg(sprintf("  weight matrix %d x %d in %s", nrow(wm), ncol(wm), took(t0)))

t0 <- Sys.time()
link_list <- as.data.table(getLinkList(wm))
log_msg(sprintf("  %d edges extracted in %s", nrow(link_list), took(t0)))

# getLinkList returns regulatoryGene/targetGene/weight already, but the column names have
# changed across GENIE3 versions ("regulatoryGene"/"fromNode"), so normalise by position
# after asserting the shape rather than trusting the names.
if (ncol(link_list) != 3) {
  stop(sprintf("unexpected getLinkList shape: %d columns (%s)",
               ncol(link_list), paste(names(link_list), collapse = ", ")))
}
setnames(link_list, c("regulatoryGene", "targetGene", "weight"))

setorder(link_list, -weight)
if (opts$top_n > 0 && nrow(link_list) > opts$top_n) {
  log_msg(sprintf("  truncating to top %d edges by weight", opts$top_n))
  link_list <- link_list[seq_len(opts$top_n)]
}

dir.create(dirname(opts$out), showWarnings = FALSE, recursive = TRUE)
fwrite(link_list, opts$out)
log_msg(sprintf("wrote %s (%d edges, %d regulators, %d targets)",
                opts$out, nrow(link_list),
                uniqueN(link_list$regulatoryGene), uniqueN(link_list$targetGene)))
log_msg("upload this as the GENIE3 network (--genie_network)")
log_msg(sprintf("total runtime %s", took(t_start)))
