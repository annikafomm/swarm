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

# example command to run this script:
# Rscript ../backend/calc_multiome_scores/calc_multiome_scores_test.R \
#   --outdir /workspaces/swarm/frontend/../backend/uploads/job_1772776453097_2061e2f5-d083-4356-b634-aceaed13aad6/dumpase1 \
#   --log /workspaces/swarm/frontend/../backend/uploads/job_1772776453097_2061e2f5-d083-4356-b634-aceaed13aad6/dumpase1/calc_scores.log \
#   --multiome_rds /workspaces/swarm/backend/uploads/job_1772776453097_2061e2f5-d083-4356-b634-aceaed13aad6/healthy_breast_preprocessed_RNA_peaks_only.rds \
#   --fragments_tsv_gz /workspaces/swarm/backend/uploads/job_1772776453097_2061e2f5-d083-4356-b634-aceaed13aad6/fake_fragments.tsv.gz \
#   --fragments_tsv_gz_tbi /workspaces/swarm/backend/uploads/job_1772776453097_2061e2f5-d083-4356-b634-aceaed13aad6/fake_fragments.tsv.gz.tbi \
#   --multiome \
#   --chromvar \
#   --footprinting \
#   --genome hg38 \
#   --differential_motif_activity

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

#' Project chromVAR motif deviation scores from dissociated cells onto spatial spots
#'
#' ChromVAR deviation scores computed per dissociated cell are projected onto
#' spatial spots via a weighted mean using the Tangram mapping matrix.
#' The resulting spot-level scores are stored as a new \code{chromvar} assay
#' inside \code{spot_obj} and active idents are set to \code{clustering_var}.
#'
#' @param object_dissociated A Seurat object with a chromVAR assay
#'   (motifs \eqn{\times} cells deviation scores).
#' @param spot_meta A data.frame with spot metadata (e.g., from spatial data).
#'   \code{seuratObj_dissociated2spatial()}. Its column names must correspond to
#'   spatial spot IDs.
#' @param M A numeric matrix or data.frame of dimensions
#'   \code{(n_cells \times n_spots)} containing Tangram mapping probabilities.
#'   Row names must match \code{colnames(object_dissociated)}.
#' @param chromvar_assay Name of the chromVAR assay in \code{object_dissociated}
#'   (default: \code{"chromvar"}).
#' @param clustering_var Column in \code{spot_obj@meta.data} that carries
#'   spatial cluster labels; set as active idents after projection
#'   (default: \code{"leiden"}).
#'
#' @return \code{spot_obj} with a new \code{chromvar} assay containing the
#'   projected deviation scores and active idents set to \code{clustering_var}.
project_chromvar_to_spots <- function(
    object_dissociated,
    spot_meta,
    M,
    chromvar_assay = "chromvar",
    clustering_var = "leiden"
) {
    # LIKE IN NOTEBOOK
    C <- GetAssayData(object_dissociated, assay = chromvar_assay, layer = "data")
    M_sp <- Matrix::Matrix(as.matrix(M), sparse = TRUE)
    C_spot <- C %*% M_sp

    fix_barcode <- function(x) sub("\\.(\\d+)$", "-\\1", x)
    # fix spot IDs consistently
    colnames(C_spot) <- fix_barcode(colnames(C_spot))

    to_dgC <- function(x) {
      if (inherits(x, "dgCMatrix")) return(x)
      Matrix::Matrix(as.matrix(x), sparse = TRUE)
    }

    # Make sure C_spot is motifs x spots
    C_spot_dgC <- to_dgC(C_spot)
    colnames(spot_meta)[1] <- "X"
    rownames(spot_meta) <- spot_meta$X

    # Create Seurat object from the MATRIX (not Assay)
    spot_obj <- CreateSeuratObject(
      counts = C_spot_dgC,
      assay = chromvar_assay
    )

    spot_obj <- SetAssayData(
      object = spot_obj,
      assay  = "chromvar",
      layer  = "data",
      new.data = GetAssayData(spot_obj, assay = "chromvar", layer = "counts")
    )
    cat("spot_obj barcodes:", head(colnames(spot_obj)), "\n")
    cat("spot_meta rownames:", head(rownames(spot_meta)), "\n")
    # Add metadata (now aligned)
    spot_obj <- AddMetaData(spot_obj, spot_meta)
    return(spot_obj)
}


#' Differential motif activity for a given identity using chromVAR
#'
#' Uses \code{Seurat::FindMarkers} on a chromVAR assay to identify motifs with
#' differential activity between \code{ident1} and \code{ident2}. Results can
#' optionally be stored in \code{object@misc$diff_motif_activity}.
#'
#' @param object A Seurat object containing a chromVAR assay.
#' @param ident1 Identity (cluster) of interest.
#' @param ident2 Optional identity to compare against. If \code{NULL}, all
#'   remaining cells are treated as \code{"rest"} when storing in \code{misc}.
#' @param save_in_misc Logical; if \code{TRUE}, saves results under
#'   \code{object@misc$diff_motif_activity} and returns the object. If
#'   \code{FALSE}, returns the differential activity table.
#' @param chromvar_assay Name of the chromVAR assay (default: \code{"chromvar"}).
#' @param only_pos Logical; whether to return only positive markers
#'   (default: \code{TRUE}).
#'
#' @return If \code{save_in_misc = TRUE}, the updated Seurat object. Otherwise,
#'   a data frame of differential motif activity statistics.
diff_motif_activity_for_ident <- function(
    object,
    ident1,
    ident2 = NULL,
    save_in_misc = TRUE,
    chromvar_assay = "chromvar",
    only_pos = TRUE
    ){
    differential.activity <- FindMarkers(
        object = object,
        ident.1 = ident1,
        ident.2 = ident2,
        assay = chromvar_assay,
        only.pos = only_pos,
        mean.fxn = rowMeans,
        fc.name = "avg_diff"
        )
    if (save_in_misc){
        if (is.null(object@misc$diff_motif_activity)){
            object@misc$diff_motif_activity <- list()
        }
        if (is.null(ident2)){
            ident2 <- "rest"
        }
        # save under ident1_vs_ident2
        object@misc$diff_motif_activity[[paste0(as.character(ident1), "_vs_", as.character(ident2))]] <- differential.activity
        return(object)
    }
    return(differential.activity)
}

run_diff_motif_activity_spatial <- function(spot_obj,
                                            ident_col = "cell_type",
                                            chromvar_assay = "chromvar",
                                            min_n = 3,
                                            only_pos = FALSE,
                                            logfc.threshold = 0.1,
                                            save_in_misc = TRUE) {
  if (!ident_col %in% colnames(spot_obj@meta.data)) {
    stop(sprintf("ident_col='%s' not found in spot_obj@meta.data.", ident_col))
  }
  Idents(spot_obj) <- spot_obj[[ident_col, drop = TRUE]]

  idents <- levels(Idents(spot_obj))
  for (ident in idents) {
    n1 <- sum(Idents(spot_obj) == ident)
    if (n1 < min_n) next
    spot_obj <- diff_motif_activity_for_ident(
      object = spot_obj,
      ident1 = ident,
      ident2 = NULL,
      save_in_misc = save_in_misc,
      chromvar_assay = chromvar_assay,
      only_pos = only_pos
      #logfc.threshold = logfc.threshold
    )
  }
  return(spot_obj)
}


spatialObj_add_fragments_dissociated2spatial <- function(
    object_dissociated,
    M, # tangram map eg M = read.csv("../swarm/backend/uploads/job_1768478211396_76193c44-3781-4b8b-b9ef-fc7702b4a617/plasmidpoop/tangram_map.csv",row.names=1), adata_map.X
    spot_obj,
    assay = "peaks",
    slot = "counts"
){
    C <- GetAssayData(object_dissociated, assay = assay, slot = slot)
    # force to 2D numeric (and to sparse for speed)
    M_sp <- Matrix::Matrix(as.matrix(M), sparse = TRUE)

    C_spot <- C %*% M_sp
    fix_barcode <- function(x) sub("\\.(\\d+)$", "-\\1", x)
    colnames(C_spot) <- fix_barcode(colnames(C_spot))

    ranges <- granges(object_dissociated[[assay]])  # peak GRanges from your multiome object

    chrom_spot <- CreateChromatinAssay(
        counts = C_spot,
        ranges = ranges,
        genome = genome(object_dissociated[[assay]]), # or "hg38"
        motifs = Motifs(object_dissociated[[assay]]),
        annotation = Annotation(object_dissociated[[assay]]),
        #fragments = Fragments(object_dissociated[[assay]])
        # Fragments is per cell, so we cannot transfer it here/ would have to aggregate with M
    )

    # spot_obj <- CreateSeuratObject(counts = chrom_spot, assay = assay)
    # Add new assay to existing spot_obj
    spot_obj[[assay]] <- chrom_spot

    DefaultAssay(spot_obj) <- assay

    spot_obj <- AddMetaData(spot_obj, spot_obj@meta.data)
    fix_barcode <- function(x) sub("\\.(\\d+)$", "-\\1", x)
    return(spot_obj)

}

plot_footprint_for_motif <- function(
  motif_name,
  M=NULL,
  M_path=NULL, # eg adata_map.X.csv
  spot_obj=NULL,
  spot_obj_path=NULL,
  object_dissociated = NULL,
  object_dissociated_path = NULL,
  assay = "peaks",
  clustering_var = "cell_type",
  object_dissociated_out_path = "dissociated_obj_out.rds",
  spot_obj_out_path = "spot_obj_out.rds",
  plot_out_path = "test.pdf",
  overwrite = TRUE,
  save_objects = FALSE
) {

  # Safety checks for overwrite behavior
  if (!overwrite) {
    if (file.exists(object_dissociated_out_path)) {
      stop("Refusing to overwrite existing file: ", object_dissociated_out_path)
    }
    if (file.exists(spot_obj_out_path)) {
      stop("Refusing to overwrite existing file: ", spot_obj_out_path)
    }
  }
  if (is.null(object_dissociated)) {
    if (is.null(object_dissociated_path)) {
      stop("Either object_dissociated or object_dissociated_path must be provided.")
    }
    object_dissociated <- readRDS(object_dissociated_path)
  }
  if (is.null(M)) {
    if (is.null(M_path)) {
      stop("Either M or M_path must be provided.")
    }
    M <- read.csv(M_path, row.names = 1, check.names = FALSE)
  }
  if (is.null(spot_obj)) {
    if (is.null(spot_obj_path)) {
      stop("Either spot_obj or spot_obj_path must be provided.")
    }
    spot_obj <- readRDS(spot_obj_path)
  }

  key <- sprintf("footprint_%s", paste(motif_name))
  print("class of peaks assay in original object:")
  print(class(object_dissociated[["peaks"]]))

  # Compute footprint in dissociated object
  object_dissociated <- Signac::Footprint(
    object = object_dissociated,
    assay = assay,
    motif.name = motif_name,
    genome = BSgenome.Hsapiens.UCSC.hg38,
    in.peaks = TRUE,
    key = key
  )

  # Project dissociated footprint to spots
  fp_spots_all <- footprints_dissociated2spatial(
    object_dissociated = object_dissociated,
    footprint_key = key,
    M = M,
    assay = assay
  )
  fix_barcode <- function(x) sub("\\.(\\d+)$", "-\\1", x)
  rownames(fp_spots_all) <- fix_barcode(rownames(fp_spots_all))
  # Store in spot object
  spot_obj[[assay]]@positionEnrichment[[key]] <- fp_spots_all
  Seurat::Idents(spot_obj) <- clustering_var

  if (save_objects) {
    # Persist updated objects to disk
    saveRDS(object_dissociated, file = object_dissociated_out_path)
    saveRDS(spot_obj, file = spot_obj_out_path)
  }

  # Plot
  p <- Signac::PlotFootprint(
    object = spot_obj,
    features = key,
    group.by = clustering_var,
    show.expected = TRUE,
    normalization = "subtract"
  )
  if (!is.null(plot_out_path)) {
    ggsave(filename = plot_out_path, plot = p, width = 6, height = 4)
  }else{
    plot_out_path <- file.path(dirname(spot_obj_out_path),paste0("footprint_", motif_name, ".pdf"))
    ggsave(filename = plot_out_path, plot = p, width = 6, height = 4)
  }
  return(list(
    plot = p,
    object_dissociated = object_dissociated,
    spot_obj = spot_obj,
    plot_out_path = plot_out_path,
    footprint_key = key
  ))

}


# top_motifs_per_comparison <- function(spot_obj, n = 5) {
#   comparisons <- names(spot_obj@misc$diff_motif_activity)
#   top_motifs <- lapply(comparisons, function(comp) {
#     res <- spot_obj@misc$diff_motif_activity[[comp]]
#     if (is.null(res) || nrow(res) == 0) return(res)

#     res <- as.data.frame(res)
#     res$motif <- rownames(res)

#     score_col <- c("avg_diff", "avg_log2FC", "avg_logFC")[
#       c("avg_diff", "avg_log2FC", "avg_logFC") %in% colnames(res)
#     ][1]

#     if (is.na(score_col)) {
#       stop("No ranking column found for comparison: ", comp)
#     }

#     ord <- order(res[[score_col]], decreasing = TRUE, na.last = NA)
#     head(res[ord, , drop = FALSE], n = n)
#   })
#   names(top_motifs) <- comparisons
#   return(top_motifs)
# }
top_motifs_per_comparison <- function(spot_obj, n = 5) {
  comparisons <- names(spot_obj@misc$diff_motif_activity)

  top_motifs <- lapply(comparisons, function(comp) {
    res <- spot_obj@misc$diff_motif_activity[[comp]]
    if (is.null(res) || nrow(res) == 0) return(NULL)

    res <- as.data.frame(res)
    res$motif <- rownames(res)

    score_col <- c("avg_diff", "avg_log2FC", "avg_logFC")[
      c("avg_diff", "avg_log2FC", "avg_logFC") %in% colnames(res)
    ][1]

    if (is.na(score_col)) {
      stop("No ranking column found for comparison: ", comp)
    }

    ord <- order(res[[score_col]], decreasing = TRUE, na.last = NA)
    head(res[ord, , drop = FALSE], n = n)
  })

  names(top_motifs) <- comparisons
  return(top_motifs)
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
    library(BiocParallel)
    register(SerialParam())
    object <- RunChromVAR(object = object, genome = genome, assay = "peaks")

    # save chromVAR scores
    chromvar_scores <- as.data.frame(t(GetAssayData(object = object, assay = "chromvar", slot = "data")))
    outdir_multiome <- file.path(outdir, "multiome")
    if (!dir.exists(outdir_multiome)) {
            dir.create(outdir_multiome, recursive = TRUE, showWarnings = FALSE)
          }
    write.csv(chromvar_scores, file = file.path(outdir_multiome, "chromvar_scores.csv"))

    if (args$differential_motif_activity) {
      log_message("Differential motif activity on spatial spots (projected chromVAR)...", logfile, 2)

      M_path <- file.path(args$outdir, "adata_map.X.csv")
      meta_path <- file.path(args$outdir, "adata_map.var.csv")

      if (!file.exists(M_path)) stop("Tangram map not found: ", M_path)
      if (!file.exists(meta_path)) stop("Spot metadata not found: ", meta_path)

      M <- read.csv(M_path, row.names = 1, check.names = FALSE)
      spot_meta <- read.csv(meta_path, check.names = FALSE)

      spot_obj <- project_chromvar_to_spots(
        object_dissociated = object,
        M = M,
        spot_meta = spot_meta,
        chromvar_assay = "chromvar",
        clustering_var = "leiden"
      )

      spot_obj <- run_diff_motif_activity_spatial(
        spot_obj = spot_obj,
        ident_col = args$spot_groupby,
        chromvar_assay = "chromvar",
        min_n = args$diff_motif_min_spots,
        only_pos = args$diff_motif_only_pos,
        logfc.threshold = args$diff_motif_logfc_threshold,
        save_in_misc = TRUE
      )

      saveRDS(spot_obj, file.path(outdir, "spot_obj_chromvar.rds"))

      top_tbls <- top_motifs_per_comparison(spot_obj, n = 3)


      multiome_dir <- file.path(outdir, "multiome")
      if (!dir.exists(multiome_dir)) {
        dir.create(multiome_dir, recursive = TRUE)
      }

      for (comp in names(top_tbls)) {
        tbl <- top_tbls[[comp]]
        if (is.null(tbl) || nrow(tbl) == 0) next

        safe_comp <- gsub("[^A-Za-z0-9._-]", "_", comp)
        out_file <- file.path(
          multiome_dir,
          paste0("diff_motif_activity_top_motifs_", safe_comp, ".csv")
        )

        write.csv(tbl, out_file, row.names = FALSE)
      }
      #write.csv(top_tbls, file.path(outdir, "diff_motif_activity_top_motifs.csv"), row.names = FALSE)
    }

  }
  # Enable footprinting for given motifs or da motifs
  if (args$footprinting) {
    log_message("Footprinting: compute on dissociated, project to spatial spots, plot...", logfile, 2)

    # --- Inputs: Tangram map + spot metadata ---
    if (is.null(M)) {
      M_path <- file.path(args$outdir, "adata_map.X.csv")
      if (!file.exists(M_path)) stop("Tangram map not found: ", M_path)
      M <- read.csv(M_path, row.names = 1, check.names = FALSE)
    }
    if (is.null(spot_meta)) {
      meta_path <- file.path(args$outdir, "adata_map.var.csv")
      if (!file.exists(meta_path)) stop("Spot metadata not found: ", meta_path)
      spot_meta <- read.csv(meta_path, check.names = FALSE)
    }


    # --- Fragment file (required by Signac Footprint) ---
    if (is.null(args$fragments_tsv_gz) || !nzchar(args$fragments_tsv_gz)) {
      stop("Footprinting requires --fragments_tsv_gz (path to fragments.tsv.gz).")
    }
    if (!file.exists(args$fragments_tsv_gz)) {
      stop("fragments_tsv_gz not found: ", args$fragments_tsv_gz)
    }
    file_info <- file.info(args$fragments_tsv_gz)
    if (file_info$size == 0) {
        cat("Fragment file empty, using default path\n")
        args$fragments_tsv_gz <- "../backend/data/multiome_data/GSM7821196_KTBpool6_Ashkenazi_jew_atac_fragments.tsv.gz"
    }

    # Attach fragments to dissociated peaks assay (needs to match the cells in object)
    frag_obj <- CreateFragmentObject(path = args$fragments_tsv_gz, cells = colnames(object))
    object[["peaks"]]@fragments <- list(frag_obj)

    if (is.null(spot_obj)) {
      # TODO: enable Footprints w/o dma
    }
    spot_obj <- spatialObj_add_fragments_dissociated2spatial(
        object_dissociated= object,
        M=M, # tangram map eg M = read.csv("../swarm/backend/uploads/job_1768478211396_76193c44-3781-4b8b-b9ef-fc7702b4a617/plasmidpoop/tangram_map.csv",row.names=1), adata_map.X
        spot_obj=spot_obj,
        assay = "peaks",
        slot = "counts"
    )

    # --- Which motifs to footprint? ---
    motifs <- NULL
    if (!is.null(args$footprinting_motifs) && nzchar(args$footprinting_motifs)) {
      # list of motifs given as argument, comma separated
      motifs <- trimws(unlist(strsplit(args$footprinting_motifs, ",")))
    } else if (!is.null(spot_obj@misc$diff_motif_activity)) {
      # dma results available, use top motifs from there
      top_tbls <- top_motifs_per_comparison(spot_obj, n = args$footprinting_top_n)
      motifs <- c()

      for (comp in names(top_tbls)) {
          res <- top_tbls[[comp]]
          if (!is.null(res) && nrow(res) > 0 && "motif" %in% colnames(res)) {
              motifs <- c(motifs, as.character(res$motif))
          }
      }

      motifs <- unique(motifs)
      motifs <- motifs[!is.na(motifs) & nzchar(motifs)]

    } else {
      # no motifs specified and no dma results, fallback to hardcoded list (for testing)
      motifs <- c("MA0084.2")
    }

    motifs <- motifs[!is.na(motifs) & nzchar(motifs)]
    motifs <- unique(motifs)

    if (length(motifs) == 0) stop("No motifs selected for footprinting.")

    log_message(paste0("Footprinting motifs (n=", length(motifs), "): ", paste(motifs, collapse = ", ")), logfile, 2)

    # --- Run footprints and plot ---
    for (motif in motifs[0:3]) {
      plot_out_path <- file.path(outdir, paste0("footprint_", motif, ".pdf"))
      res <- plot_footprint_for_motif(
        motif_name = motif,
        M = M,
        spot_obj = spot_obj,
        object_dissociated = object,
        assay = "peaks",
        clustering_var = "cell_type",
        plot_out_path = plot_out_path,
        overwrite = TRUE
      )
      object <- res$object_dissociated
      spot_obj <- res$spot_obj
    }

    saveRDS(spot_obj, file.path(outdir, "spot_obj_footprints.rds"))
    saveRDS(object, file.path(outdir, "dissociated_obj_footprints.rds"))
    log_message("Footprinting done. Wrote spot_obj_footprints.rds and dissociated_obj_footprints.rds and footprint_*.pdf files.", logfile, 2)
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



    # tangram / spatial inputs (for spatial projected chromVAR + footprinting)
    make_option("--tangram_map_csv", type="character", default="", help="Path to Tangram mapping matrix CSV (cells x spots). Default: <outdir>/adata_map.X.csv", metavar="file"),
    make_option("--spot_meta_csv", type="character", default="", help="Path to spot metadata CSV. Default: <outdir>/adata_map.var.csv", metavar="file"),
    make_option("--spot_meta_id_col", type="character", default="X", help="Column in spot_meta_csv containing spot barcodes if rownames are not set (often 'X').", metavar="col"),
    make_option("--spot_groupby", type="character", default="cell_type", help="Metadata column used as identities for spatial differential motif activity and footprinting.", metavar="col"),

    # differential motif activity params (spatial)
    make_option("--diff_motif_min_spots", type="integer", default=3, help="Minimum number of spots per identity to run differential motif activity.", metavar="int"),
    make_option("--diff_motif_logfc_threshold", type="double", default=0.1, help="FindMarkers logfc.threshold for differential motif activity (spatial).", metavar="float"),
    make_option("--diff_motif_only_pos", action="store_true", default=FALSE, help="If set, only report positive markers in differential motif activity (spatial)."),

    # footprinting params
    make_option("--footprinting_motifs", type="character", default="", help="Comma-separated motif IDs to footprint. If empty, uses top motifs from diff_motif_activity_top_motifs.csv / spot_obj@misc when available.", metavar="str"),
    make_option("--footprinting_top_n", type="integer", default=3, help="Number of top motifs per comparison to use when --footprinting_motifs is empty.", metavar="int"),
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
