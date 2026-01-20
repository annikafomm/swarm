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

# source("../backend/calc_R_scores/utils.R")
source("../backend/calc_multiome_scores/calc_global_motif_analysis.R")

#  Rscript calc_multiome_scores.R --outdir path/to/out_dir -log path/to/log_file -multiome -chromvar --multiome_rds path/to/multiome_rds path/to/multiome_rds_file
# eg:
# Rscript ../backend/calc_multiome_scores/calc_multiome_scores.R --outdir /workspaces/swarm/frontend/../backend/uploads/job_1764602732819_merit/dumpase1 --log /workspaces/swarm/frontend/../backend/uploads/job_1764602732819_merit/dumpase1/calc_scores.log --multiome_rds /workspaces/swarm/backend/uploads/job_1764602732819_merit/healthy_breast_preprocessed_RNA_peaks_only.rds --multiome --chromvar --genome hg38
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

  print("Running global motif analysis...")

  print("checking frag path")
  fr <- Fragments(object[["peaks"]])
  # what Signac thinks the fragments file is
  old_path <- fr[[1]]@path
  print(old_path)

  outdir <- args$outdir
  print("args$chromvar:")
  print(args$chromvar)


  if (args$chromvar | args$differential_motif_activity | args$motif_enrichment | args$footprinting ) {
    if (args$genome == "hg38") {
      suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg38))
      genome <- BSgenome.Hsapiens.UCSC.hg38
    } else if (args$genome %in% c("hg19", "hg37")) {
      suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg19))
      genome <- BSgenome.Hsapiens.UCSC.hg19
    } else {
      stop("Unknown genome: ", args$genome, " (use hg38 or hg19/hg37)")
    }
    log_message("Adding Motifs ...", logfile, 2)
    object <- add_jaspar2024_motifs(object, out_path=outdir, genome=genome)
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



  if (args$chromvar) {
    # run chromVAR
    # object[["chromvar"]]@data      # motif deviation scores (motifs × cells)
    # object[["chromvar"]]@meta.features  # motif-level metadata, if present
    log_message("Running chromVAR...", logfile, 2)
    object <- RunChromVAR(object = object, genome = genome, assay = "peaks")

    # save chromVAR scores
    chromvar_scores <- as.data.frame(t(GetAssayData(object = object, assay = "chromvar", slot = "data")))
    outdir_multiome <- file.path(outdir, "multiome")
    if (!dir.exists(outdir_multiome)) {
            dir.create(outdir_multiome, recursive = TRUE, showWarnings = FALSE)
          }
    write.csv(chromvar_scores, file = file.path(outdir_multiome, "chromvar_scores.csv"))

    if (args$differential_motif_activity) {
      # find differential motif activity for each cluster
      print("Finding differential motif activity...")
      idents <- levels(object)
      for (ident in idents){
          # save in object@misc$diff_motif_activity[[paste0(as.character(ident1), "_vs_", as.character(ident2))]]
          object <- diff_motif_activity_for_ident(object, ident1 = ident, save_in_misc = TRUE)
      }
    }

  }

  if (args$footprinting) {
    log_message("Running footprinting analysis...", logfile, 2)

    log_message("setting fragment path...", logfile, 2)
    frag_obj <- CreateFragmentObject(
      path  = args$fragments_tsv_gz,
      cells = colnames(object)   # restrict to cells in this object
    )
    object[["peaks"]]@fragments <- list(frag_obj)
    log_message(Fragments(object[["peaks"]])[[1]]@path, logfile, 2)

    # compute Tn5 insertion bias for footprinting, is saved in object[["peaks"]]@bias
    log_message("Computing Tn5 insertion bias...", logfile, 2)
    object <- compute_Tn5_insertion_Bias(object, genome = genome, assay = "peaks")
    # create spatial seurat object
    M_path <- file.path(args$outdir, "adata_map.X.csv")
    M <- read.csv(M_path, row.names = 1)
    spot_meta <- read.csv(file.path(args$outdir, "adata_map.var.csv"))
    spot_obj <- seuratObj_dissociated2spatial(
        object_dissociated= object,
        M=M,
        spot_meta=spot_meta,
        assay = "peaks",
        slot = "counts")
    spot_obj_out_path <- file.path(args$outdir, "spot_obj.rds")
    saveRDS(spot_obj, spot_obj_out_path)

    # plot footprint for top motifs
    log_message("Plotting footprints for top motifs...", logfile, 2)
    object_dissociated_path <- args$multiome_rds
    motif_name <- "MA0084.2"
    plot_footprint_for_motif(
        M_path = M_path, # eg adata_map.X.csv
        spot_obj_path = spot_obj_out_path,
        motif_name = motif_name,
        object_dissociated = multiome_data,
        object_dissociated_path = object_dissociated_path, # wont be used since object_dissociated is provided
        assay = "peaks",
        clustering_var = "leiden",
        object_dissociated_out_path = object_dissociated_path,
        spot_obj_out_path = spot_obj_out_path,
        plot_out_path = NULL,
        overwrite = TRUE)
  }
  return(object)
}




main <- function() {
  option_list <- list(
    # dir with files
    make_option("--multiome_rds", type="character", help="Path to multiome RDS file", metavar="file"),
    make_option("--outdir", type="character", help="Dir containing output files"),
    make_option("--log", type="character", help="Path to the log file", metavar="file"),
    make_option("--fragments_tsv_gz", type="character", help="Path to fragments TSV.GZ file for footprinting", metavar="file"),
    make_option("--fragments_tsv_gz_tbi", type="character", help="Path to fragments TSV.GZ.TBI file for footprinting", metavar="file"),


    # score flags
    make_option("--chromvar", action="store_true", default=FALSE, help="Calculate chromVAR score"),
    make_option("--differential_motif_activity", action="store_true", default=FALSE, help="Calculate differential motif activity"),
    make_option("--motif_enrichment", action="store_true", default=FALSE, help="Calculate motif enrichment for idents"),
    make_option("--footprinting", action="store_true", default=FALSE, help="Calculate footprinting Tn5 insertion bias"),

    # score flags
    make_option("--multiome", action="store_true", default=FALSE, help="Calculate for multiome"),

    # motif mapping flag
    make_option("--genome", type = "character", default = "hg38", help = "reference genome: hg38 or hg19 [default: %default]")
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
  print("logfile")
  print(logfile)

  log_message(paste0("R Multiome Pipeline started at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), logfile)
  print("multiome params:")
  print(args)

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
