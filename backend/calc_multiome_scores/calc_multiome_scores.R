suppressPackageStartupMessages({
  # suppressWarnings(library(Matrix))
  # suppressWarnings(library(data.table))
  # suppressWarnings(library(viper))
  # suppressWarnings(library(jsonlite))
  # suppressWarnings(library(optparse))
  # suppressWarnings(library(SPONGE))
  # suppressWarnings(library(dplyr))
  # suppressWarnings(library(tibble))
  # suppressWarnings(library(GSVA))
  # suppressWarnings(library(AUCell))
  suppressWarnings(library(Signac))
  suppressWarnings(library(Seurat))
  suppressWarnings(library(JASPAR2024))
  suppressWarnings(library(TFBSTools))
  suppressWarnings(library(BSgenome.Hsapiens.UCSC.hg38))
  suppressWarnings(library(patchwork))
  suppressWarnings(library(ggplot2))
})

# source("../backend/calc_R_scores/utils.R")
source("../backend/calc_multiome_scores/calc_global_motif_analysis.R")

#  Rscript calc_scores.R   --dir ../datasets_scores/GSM6592049_M2_prepro   --dir ../datasets_scores/GSM6592049_M2_scores   --sponge_network ../networks/SPONGE/breast_invasive_carcinoma/breast_invasive_carcinoma_interactionNetwork.csv   --sponge_analysis ../networks/SPONGE/breast_invasive_carcinoma/breast_invasive_carcinoma_networkAnalysis.csv --genie_network ../networks/GENIE3/BRCA/genie3_BRCA_mrn.top_100k.csv --ensembl_col ensemble_id  --aucell --gsva --ssgsea --viper
#  Rscript calc_multiome_scores.R --outdir path/to/out_dir -log path/to/log_file -multiome -chromvar --multiome_rds path/to/multiome_rds path/to/multiome_rds_file

# ! parallel::detectCores()

# logging helper
log_message <- function(msg, logfile, indent = 0) {
  prefix <- paste(rep(" ", indent), collapse = "")
  line <- paste0(prefix, msg)   # <-- no "\n" here
  cat(line, "\n", sep = "")     # console
  write(line, file = logfile, append = TRUE)  # file (no extra newline)
}

format_runtime <- function(t0) {
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  mins <- floor(elapsed / 60)
  secs <- round(elapsed %% 60)
  sprintf("%d min %d sec", mins, secs)
}

global_motif_analysis <- function(object, args, logfile) {

  outdir <- args$outdir

  if (args$chromvar | args$differential_motif_activity | args$motif_enrichment) {
    log_message("Adding Motifs ...", logfile, 2)
    object <- add_jaspar2024_motifs(object, genome=genome)
  }

  if (args$motif_enrichment){
    idents <- levels(object)
    log_message("Finding differentially accessible peaks...", logfile, 2)
    for (ident in idents){
        # saved in object@misc$diff_accessible_peaks[[paste0(as.character(ident1), "_vs_", as.character(ident2))]]
        object <- diff_peaks_for_ident(object, ident1 = ident, save_in_misc = TRUE)
    }
    # find enriched motifs
    # saved in object@misc$enriched_motifs[[comparison]] with comparison in names(object@misc$diff_accessible_peaks)
    log_message("Finding enriched motifs...", logfile, 2)
    object <- find_enriched_motifs(object, save_in_misc = TRUE)
  }

  if (args$footprinting | args$chromvar) {
    if (opt$genome == "hg38") {
      suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg38))
      genome <- BSgenome.Hsapiens.UCSC.hg38
    } else if (opt$genome %in% c("hg19", "hg37")) {
      suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg19))
      genome <- BSgenome.Hsapiens.UCSC.hg19
    } else {
      stop("Unknown genome: ", opt$genome, " (use hg38 or hg19/hg37)")
    }
  }

  if (args$chromvar) {
    # run chromVAR
    # object[["chromvar"]]@data      # motif deviation scores (motifs × cells)
    # object[["chromvar"]]@meta.features  # motif-level metadata, if present
    log_message("Running chromVAR...", logfile, 2)
    object <- RunChromVAR(object = object, genome = genome, assay = "peaks")

    # save chromVAR scores
    chromvar_scores <- as.data.frame(t(GetAssayData(object = object, assay = "chromvar", slot = "data")))
    write.csv(chromvar_scores, file = file.path(outdir, "chromvar_scores.csv"))

    if (args$differential_motif_activity) {
      # find differential motif activity for each cluster
      print("Finding differential motif activity...")
      idents <- levels(object)
      for (ident in idents){
          # save in object@misc$diff_motif_activity[[paste0(as.character(ident1), "_vs_", as.character(ident2))]]
          object <- diff_motif_activity_for_ident(object, ident1 = ident, save_in_misc = TRUE)
      }
    }

    if (args$footprinting) {
      # compute Tn5 insertion bias for footprinting, is saved in object[["peaks"]]@bias
      log_message("Computing Tn5 insertion bias...", logfile, 2)
      object <- compute_Tn5_insertion_Bias(object, genome = genome, assay = "peaks")
    }
  }
  return(object)
}


main <- function() {
  option_list <- list(
    # dir with files
    make_option("--multiome_rds", type="character", help="Path to multiome RDS file", metavar="file"),
    make_option("--outdir", type="character", help="Dir containing output files"),
    make_option("--log", type="character", help="Path to the log file", metavar="file"),

    # score flags
    make_option("-chromvar", action="store_true", default=FALSE, help="Calculate chromVAR score"),
    make_option("-differential_motif_activity", action="store_true", default=FALSE, help="Calculate differential motif activity"),
    make_option("-motif_enrichment", action="store_true", default=FALSE, help="Calculate motif enrichment for idents"),
    make_option("-footprinting", action="store_true", default=FALSE, help="Calculate footprinting Tn5 insertion bias"),
    # score flags
    make_option("-multiome", action="store_true", default=FALSE, help="Calculate for multiome"),

    # motif mapping flag
    make_option("--genome", type = "character", default = "hg38", help = "reference genome: hg38 or hg19 [default: %default]"),
  )

  parser <- OptionParser(option_list=option_list)
  args <- parse_args(parser)

  # required flags
  if (is.null(args$log)) {
    stop("Error: --log is required.")
  }
  if (is.null(args$outdir)) {
    stop("Error: --outdir is required.")
  }
  if (is.null(args$multiome_rds)) {
    stop("Error: --multiome_rds is required.")
  }

  # Prepare log file
  logfile <- args$log
  log_message(paste0("R Multiome Pipeline started at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), logfile)

  # Load multiome RDS
  log_message("Loading multiome RDS file ...", logfile, 2)
  t0 <- Sys.time()
  multiome_data <- readRDS(args$multiome_rds)
  log_message(sprintf("Multiome RDS file loaded in %s", format_runtime(t0)), logfile, 4)


  # Global Motif analysis steps
  multiome_data <- global_motif_analysis(multiome_data, args, logfile)

  out_path <- file.path(args$outdir, "global_motif_analysis.rds")
  saveRDS(multiome_data, out_path)

  log_message(paste0("R Multiome Pipeline finished at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n"), logfile)
}



# only run main() if the script is executed directly, not sourced
if (identical(environment(), globalenv()) && !length(sys.frames())) {
  main()
}
