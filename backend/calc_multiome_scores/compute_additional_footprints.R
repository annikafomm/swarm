suppressPackageStartupMessages({
  suppressWarnings(library(Matrix))
  suppressWarnings(library(data.table))
  suppressWarnings(library(viper))
  suppressWarnings(library(jsonlite))
  suppressWarnings(library(optparse))
  suppressWarnings(library(SPONGE))
  suppressWarnings(library(dplyr))
  suppressWarnings(library(tibble))
  suppressWarnings(library(GSVA))
  suppressWarnings(library(AUCell))
  suppressWarnings(library(Signac))
  suppressWarnings(library(Seurat))
  suppressWarnings(library(JASPAR2024))
  suppressWarnings(library(TFBSTools))
  suppressWarnings(library(BSgenome.Hsapiens.UCSC.hg38))
  suppressWarnings(library(patchwork))
  suppressWarnings(library(ggplot2))
})

# Determine the directory of this script so sibling files can be sourced
# regardless of the working directory the script is invoked from.
# sys.frame(1)$ofile works when source()'d; for Rscript use commandArgs.
script_dir <- tryCatch({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    dirname(normalizePath(sub("^--file=", "", file_arg[1])))
  } else {
    dirname(normalizePath(sys.frame(1)$ofile))
  }
}, error = function(e) getwd())
#source(file.path(script_dir, "calc_global_motif_analysis.R"))
source(file.path(script_dir, "calc_multiome_scores_test.R"))

# ---------------------------------------------------------------------------
# Usage example:
#   Rscript compute_additional_footprints.R \
#     --outdir /path/to/job/subdir \
#     --motif  MA1638.2 \
#     --cluster_by cell_type
# ---------------------------------------------------------------------------
option_list <- list(
  make_option("--outdir",          type = "character", help = "Job output directory (contains RDS files and adata_map.X.csv)"),
  make_option("--motif",           type = "character", default = NULL, help = "JASPAR motif ID to compute footprint for (e.g. MA1638.2)"),
  make_option("--motifs",          type = "character", default = NULL,
              help = "Comma-separated list of JASPAR motif IDs to compute footprints for (e.g. MA1638.2,MA0006.2). Overrides --motif if provided."),
  make_option("--cluster_by",      type = "character", default = "cell_type",
              help = "Metadata column to group cells by in the footprint plot [default: cell_type]"),
  make_option("--cluster_by_mult", type = "character", default = NULL,
              help = "Comma-separated list of cluster types to cluster by, overrides --cluster_by (e.g. cell_type,leiden_cluster)"),
  make_option("--save_rds",        action = "store_true", default = FALSE,
              help = "Whether to write updated RDS files with computed footprints (overwrites existing RDS files in outdir) [default: FALSE]")
)

args <- parse_args(OptionParser(option_list = option_list))

if (is.null(args$outdir)) {
  stop("--outdir is a required argument.")
}
if (is.null(args$motifs) && is.null(args$motif)) {
  stop("Either --motifs or --motif must be provided.")
}

outdir <- args$outdir
motif  <- args$motif
cluster_by <- args$cluster_by
if (!is.null(args$motifs)) {
  motifs <- unlist(strsplit(args$motifs, ","))
} else {
  motifs <- motif
}
if (!is.null(args$cluster_by_mult)) {
  cluster_by_mult <- unlist(strsplit(args$cluster_by_mult, ","))
} else {
  cluster_by_mult <- cluster_by
}


message(sprintf("[compute_additional_footprints] outdir=%s  motif=%s  cluster_by=%s",
            outdir, motif, cluster_by))

message("Checking for required files...")
dissociated_rds <- file.path(outdir, "dissociated_obj_footprints.rds")
M_path <- file.path(outdir, "adata_map.X.csv")
spot_rds       <- file.path(outdir, "spot_obj_footprints.rds")


if (!file.exists(spot_rds))        stop("spot_obj_footprints.rds not found in: ", outdir)
if (!file.exists(dissociated_rds)) stop("dissociated_obj_footprints.rds not found in: ", outdir)
if (!file.exists(M_path))          stop("Tangram map not found in: ", outdir)

message("Reading required files...")
message(paste0("spot_rds: ", spot_rds))
message(paste0("dissociated_rds: ", dissociated_rds))
message(paste0("adata_map.X.csv: ", M_path))
message("reading dissociated object...")
object   <- readRDS(dissociated_rds)
message(paste0("M_path: ", M_path))
M <- read.csv(M_path, row.names = 1, check.names = FALSE)
message("reading spot object...")
spot_obj <- readRDS(spot_rds)


available_cts <- unique(spot_obj@meta.data$cell_type)
message(paste0("Available cell types in spot_obj: ", paste(available_cts, collapse = ", ")))
# handle clusters: is_celltype: needs to be added in obs columns
for (c in cluster_by_mult) {
  if (!(c %in% colnames(spot_obj@meta.data))) {
    if (startsWith(c, "is_")) {
      ct <- sub("^is_", "", c)
      message(paste0("Processing cluster column: ", c, " with cell type: ", ct))
      if (ct %in% available_cts) {
        message(paste0("Adding column ", c, " to metadata"))
        spot_obj@meta.data[[c]] <- ifelse(spot_obj@meta.data$cell_type == ct, ct, "other cells")
      }
    }
  }
}
message(paste0("cluster_by ", cluster_by))
cluster_by_mult <- Filter(function(cl) cl %in% colnames(spot_obj@meta.data), unique(c(cluster_by_mult, cluster_by)))
message(paste0("cluster_by_mult: ", paste(cluster_by_mult, collapse = ", ")))
message(paste0("Available metadata columns in spot_obj: ", paste(colnames(spot_obj@meta.data), collapse = ", ")))


if (length(cluster_by_mult) == 0) {
  stop("None of the specified cluster types are valid metadata columns in spot_obj. Please check your --cluster_by_mult argument and ensure the specified columns exist in the spot object metadata.")
}


message("Computing footprint...")

# Optional override from environment, e.g. FOOTPRINT_WORKERS=4
n_workers <- suppressWarnings(as.integer(Sys.getenv("FOOTPRINT_WORKERS", NA)))

if (is.na(n_workers) || n_workers < 1L) {
  n_workers <- parallel::detectCores(logical = FALSE)
  if (is.na(n_workers)) n_workers <- 1L
  n_workers <- max(1L, min(length(motifs), n_workers - 1L))
} else {
  n_workers <- min(length(motifs), n_workers)
}

message(sprintf(
  "Parallelising over motifs with %d worker(s) for %d motif(s)",
  n_workers, length(motifs)
))

results <- parallel::mclapply(
  X = motifs,
  FUN = function(m) {
    tryCatch({
      for (c in cluster_by_mult) {
        message(sprintf("Processing motif: %s with cluster_by: %s", m, c))
        plot_out_path <- file.path(outdir, paste0("footprint_", m, "_", c, ".pdf"))

        plot_footprint_for_motif(
          motif_name            = m,
          M                     = M,
          spot_obj              = spot_obj,
          object_dissociated    = object,
          assay                 = "peaks",
          clustering_var        = c,
          plot_out_path         = plot_out_path,
          overwrite             = TRUE
        )
      }

      list(ok = TRUE, motif = m)
    }, error = function(e) {
      list(ok = FALSE, motif = m, error = conditionMessage(e))
    })
  },
  mc.cores = n_workers,
  mc.preschedule = FALSE
)

failed <- Filter(function(x) !isTRUE(x$ok), results)
if (length(failed) > 0L) {
  stop(
    paste(
      vapply(
        failed,
        function(x) sprintf("motif %s failed: %s", x$motif, x$error),
        character(1)
      ),
      collapse = "\n"
    )
  )
}

# message("Footprint computation complete. Saving updated objects and plot...")
# if (args$save_rds) {
#   saveRDS(object,   file.path(outdir, "dissociated_obj_footprints.rds"))
#   saveRDS(spot_obj, file.path(outdir, "spot_obj_footprints.rds"))
# }
