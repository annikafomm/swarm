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
  make_option("--outdir",     type = "character", help = "Job output directory (contains RDS files and adata_map.X.csv)"),
  make_option("--motif",      type = "character", help = "JASPAR motif ID to compute footprint for (e.g. MA1638.2)"),
  make_option("--cluster_by", type = "character", default = "cell_type",
              help = "Metadata column to group cells by in the footprint plot [default: cell_type]")
)

args <- parse_args(OptionParser(option_list = option_list))

if (is.null(args$outdir) || is.null(args$motif)) {
  stop("--outdir and --motif are required arguments.")
}

outdir <- args$outdir
motif  <- args$motif
cluster_by <- args$cluster_by


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


message("Computing footprint...")
plot_out_path <- file.path(outdir, paste0("footprint_", motif, ".pdf"))

res <- plot_footprint_for_motif(
  motif_name            = motif,
  M                     = M,
  spot_obj              = spot_obj,
  object_dissociated    = object,
  assay                 = "peaks",
  clustering_var        = cluster_by,
  plot_out_path         = plot_out_path,
  overwrite             = TRUE
)

object   <- res$object_dissociated
spot_obj <- res$spot_obj

message("Footprint computation complete. Saving updated objects and plot...")
# Persist updated objects so subsequent calls benefit from cached computations
# saveRDS(object,   file.path(outdir, "dissociated_obj_footprints_add_motifs.rds"))
# saveRDS(spot_obj, file.path(outdir, "spot_obj_footprints_add_motifs.rds"))

message(sprintf("[compute_additional_footprints] Done. Plot saved to: %s", plot_out_path))
