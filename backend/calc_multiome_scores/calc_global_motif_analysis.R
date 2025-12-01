library(Signac)
library(Seurat)
library(JASPAR2024)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(patchwork)
library(ggplot2)

# Rscript global_motif_analysis.R --in_path /nfs/.../healthy_breast_preprocessed.rds --out_path /nfs/.../healthy_breast_motif_analysis.rds
# Rscript global_motif_analysis.R --i /nfs/.../healthy_breast_preprocessed.rds --o /nfs/.../healthy_breast_motif_analysis.rds


#' Add JASPAR2024 motifs to a Seurat/Signac object
#'
#' This function adds transcription factor position frequency matrices (PFMs)
#' from the JASPAR 2024 database to an ATAC assay in a Seurat object using Signac.
#' Motifs are only added if the specified ATAC assay does not already contain motifs.
#'
#' @param object A Seurat object containing an ATAC assay.
#' @param ATAC_assay Character scalar. Name of the ATAC assay in \code{object}
#'   (default: \code{"peaks"}).
#' @param species Character scalar specifying the species used to query JASPAR
#'   (default: \code{"Homo sapiens"}).
#' @param genome A \code{BSgenome} object corresponding to the reference genome
#'   used to call peaks (default: \code{BSgenome.Hsapiens.UCSC.hg38}).
#'
#' @return The input Seurat object with motifs added to \code{object[[ATAC_assay]]}
#'   if they were not present already.
add_jaspar2024_motifs <- function(
    object,
    ATAC_assay = "peaks",
    species = "Homo sapiens",
    genome = BSgenome.Hsapiens.UCSC.hg38){
    if (is.null(object[[ATAC_assay]]@motifs)){
        jaspar <- JASPAR2024()
        sq24 <- RSQLite::dbConnect(RSQLite::SQLite(), db(jaspar))
        pfm <- TFBSTools::getMatrixSet(sq24, list(species = species, collection = "CORE", all_versions = FALSE, matrixtype = "PFM"))
        object <- AddMotifs(object, genome = genome, pfm = pfm, assay=ATAC_assay)
    }
    return(object)
}


#' Find differentially accessible peaks for a given identity
#'
#' For a given identity (cluster or group), this function computes differentially
#' accessible peaks using \code{Seurat::FindMarkers} on the ATAC assay and
#' optionally stores the results in \code{object@misc$diff_accessible_peaks}.
#'
#' @param object A Seurat object with an ATAC assay.
#' @param ident1 Identity (cluster) of interest; passed to \code{FindMarkers}
#'   as \code{ident.1}.
#' @param ident2 Optional identity to compare against. If \code{NULL}, all
#'   remaining cells are used as background and recorded as \code{"rest"} in
#'   the misc slot.
#' @param ATAC_assay Name of the ATAC assay (default: \code{"peaks"}).
#' @param test_to_use Statistical test to use in \code{FindMarkers}
#'   (default: \code{"wilcox"}).
#' @param min.pct Minimum fraction of cells in which a peak is detected,
#'   passed to \code{FindMarkers} (default: 0.05).
#' @param only_pos Logical; whether to return only positive markers
#'   (default: \code{TRUE}).
#' @param save_in_misc Logical; if \code{TRUE}, results are stored in
#'   \code{object@misc$diff_accessible_peaks} and the Seurat object is returned.
#'   If \code{FALSE}, the marker table is returned.
#' @param p_val_adj_th Adjusted p-value threshold used to filter peaks
#'   (default: 0.05).
#' @param pct.1_th Threshold on \code{pct.1} (fraction of cells in group 1)
#'   used to filter peaks (default: 0.2).
#'
#' @return If \code{save_in_misc = TRUE}, the updated Seurat object with
#'   \code{diff_accessible_peaks} stored in \code{object@misc}. Otherwise,
#'   a data frame of differentially accessible peaks.
diff_peaks_for_ident <- function(
    object,
    ident1,
    ident2 = NULL,
    ATAC_assay = "peaks",
    test_to_use = 'wilcox',
    min.pct = 0.05,
    only_pos = TRUE,
    save_in_misc = TRUE,
    p_val_adj_th = 0.05,
    pct.1_th = 0.2){

    DefaultAssay(object) <- "peaks"
    diff_accessible_peaks <- FindMarkers(
        object = object,
        ident.1 = ident1,
        ident.2 = ident2,
        only.pos = only_pos,
        test.use = test_to_use,
        min.pct = min.pct
    )
    diff_accessible_peaks <- diff_accessible_peaks[diff_accessible_peaks$p_val_adj < p_val_adj_th & diff_accessible_peaks$pct.1 > pct.1_th, ]

    if (save_in_misc){
        if (is.null(object@misc$diff_accessible_peaks)){
            object@misc$diff_accessible_peaks <- list()
        }
        if (is.null(ident2)){
            ident2 <- "rest"
        }
        # save under ident1_vs_ident2
        object@misc$diff_accessible_peaks[[paste0(as.character(ident1), "_vs_", as.character(ident2))]] <- diff_accessible_peaks
        return(object)
    }
    return(diff_accessible_peaks)
}


#' Find enriched motifs in differentially accessible peaks
#'
#' This function uses peak-level differential accessibility results stored in
#' \code{object@misc$diff_accessible_peaks} and runs \code{FindMotifs} for each
#' comparison to identify enriched transcription factor motifs.
#'
#' @param object A Seurat object with an ATAC assay and differential peak
#'   results stored in \code{object@misc$diff_accessible_peaks}.
#' @param ATAC_assay Name of the ATAC assay to use (default: \code{"peaks"}).
#' @param p_val_th P-value threshold for filtering enriched motifs
#'   (default: 0.05).
#' @param save_in_misc Logical; if \code{TRUE}, stores the list of enriched
#'   motifs in \code{object@misc$enriched_motifs} and returns the object.
#'   If \code{FALSE}, returns the list of data frames.
#'
#' @return Either the updated Seurat object (if \code{save_in_misc = TRUE}) or
#'   a named list of data frames containing enriched motifs per comparison.
find_enriched_motifs <- function(
    object,
    ATAC_assay = "peaks",
    p_val_th = 0.05,
    save_in_misc = TRUE){

    enriched_motifs_list <- list()
    diff_peaks_list <- object@misc$diff_accessible_peaks
    if (is.null(diff_peaks_list)){
        stop("No differentially accessible peaks found in object@misc$diff_accessible_peaks")
    }
    print(paste("number of diff peaks comparisons:", length(names(diff_peaks_list))))
    for (comparison in names(diff_peaks_list)){
        diff_peaks <- diff_peaks_list[[comparison]]
        peaks <- rownames(diff_peaks)
        enriched_motifs <- FindMotifs(
            object = object,
            features = peaks,
            assay = ATAC_assay
        )
        enriched_motifs <- enriched_motifs[enriched_motifs$pvalue < p_val_th, ]
        enriched_motifs_list[[comparison]] <- enriched_motifs
    }
    if (save_in_misc){
        if (is.null(object@misc$enriched_motifs)){
            object@misc$enriched_motifs <- list()
        }
        object@misc$enriched_motifs <- c(
            object@misc$enriched_motifs,
            enriched_motifs_list
        )
        return(object)
    }
    return(enriched_motifs_list)
}


#' Convert peak name strings to a GRanges object
#'
#' Converts peak names (e.g. \code{"chr1-1000-2000"} or \code{"chr1:1000-2000"}
#' depending on your format) into a \code{GRanges} object for downstream overlap
#' and genomic operations.
#'
#' @param peak_names A character vector of peak names, each encoding chromosome
#'   and coordinates. This implementation assumes a \code{"chr-start-end"}
#'   format separated by \code{"-"}.
#'
#' @return A \code{GRanges} object with one range per peak name.
peak_names_to_GRanges <- function(peak_names){

    if (length(peak_names) == 1L && grepl("'", peak_names)) {
        # split "'chr7-...''chr3-...'" into pieces
        peak_names <- strsplit(gsub("^'|'$", "", peak_names), "'\\s*'")[[1]]
    }
    parts <- do.call(rbind, strsplit(peak_names, "-"))
    gr <- GRanges(
        seqnames = parts[, 1],
        ranges   = IRanges(start = as.integer(parts[, 2]),
                        end   = as.integer(parts[, 3]))
    )
    return(gr)
}


#' Get genomic positions of a motif within differentially accessible peaks
#'
#' For a given motif and cluster comparison, this function intersects motif
#' occurrences with differentially accessible peaks stored in
#' \code{object@misc$diff_accessible_peaks}.
#'
#' @param object A Seurat object with an ATAC assay containing motif positions.
#' @param motif_name Character scalar; name of the motif in
#'   \code{object[[ATAC_assay]]$motifs@positions}.
#' @param ident1 Identity (cluster) of interest used in the differential peaks
#'   comparison.
#' @param ident2 Identity (cluster) used as background in the comparison
#'   (default: \code{"rest"}).
#' @param ATAC_assay Name of the ATAC assay containing motif information
#'   (default: \code{"peaks"}).
#'
#' @return A \code{GRanges} object of motif instances overlapping the
#'   differentially accessible peaks for the specified comparison.
get_motifs_in_da_peaks_for_cluster <- function(
    object,
    motif_name,
    ident1,
    ident2 = "rest",
    ATAC_assay = "peaks"){
        top_da_peak_gr <- peak_names_to_GRanges(
            rownames(object@misc$diff_accessible_peaks[[paste0(as.character(ident1), "_vs_", as.character(ident2))]]))

        motif_in_peaks <- subsetByOverlaps(
            object[[ATAC_assay]]$motifs@positions[[motif_name]],
            top_da_peak_gr,
            type = "any")
        return(motif_in_peaks)
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


# ' Compute Tn5 insertion bias for footprinting
# ' This function computes Tn5 insertion bias using Signac's
#' \code{InsertionBias} function and adds the bias information to the
#' Seurat object for use in footprinting analyses.
#' @param object A Seurat object with an ATAC assay.
#' @param genome A \code{BSgenome} object corresponding to the reference genome.
#' @return The updated Seurat object with Tn5 insertion bias computed.
compute_Tn5_insertion_Bias <- function(object, genome, assay = "peaks"){
    DefaultAssay(object) <- assay
    # use entire forst chromosomes as regions
    region.end <- seqlengths(x = genome)[1]
    object <- InsertionBias(
        object = object,
        genome = genome,
        region = paste0(
            names(x = region.end),
            "-1-",
            as.character(x = region.end)
        )
    )
    return(object)
}


# ' Get motif activity scores from chromVAR assay
# ' This function retrieves motif activity scores for a specified motif
#' from the chromVAR assay in a Seurat object.
#' @param object A Seurat object containing a chromVAR assay.
#' @param motif_name Character scalar; name of the motif to retrieve.
#' @param chromvar_assay Name of the chromVAR assay (default: \code{"chromvar"}).
#' @param as_numeric Logical; if \code{TRUE}, returns scores as a numeric vector
#'   (default: \code{TRUE}).
#' @return A numeric vector of motif activity scores if \code{as_numeric = TRUE},
#'   otherwise a data frame.
get_motif_activity <- function(object, motif_name, chromvar_assay = "chromvar", as_numeric = TRUE){
    DefaultAssay(object) <- chromvar_assay
    motif_activity <- FetchData(
        object = object,
        vars = motif_name,
    )
    if (as_numeric){
        motif_activity <- as.numeric(motif_activity[[motif_name]])
    }
    return(motif_activity)
}


#' Run global motif analysis pipeline
#'
#' This is a high-level convenience function that runs the complete global
#' motif analysis workflow: adding motifs, finding differentially accessible
#' peaks, computing enriched motifs, running chromVAR, and computing
#' differential motif activities. The updated object is saved to disk.
#'
#' @param in_path Path to an input RDS file containing a preprocessed Seurat
#'   multiome object.
#' @param out_path Path to an output RDS file where the updated object with
#'   motif analysis results will be saved.
#'
#' @return Invisibly returns \code{NULL}. The main side effect is writing the
#'   updated Seurat object to \code{out_path}.
main <- function(in_path, out_path, add_motifs=TRUE){
    object <- readRDS(in_path)

    # add motifs
    # saved in object[["peaks"]]@motifs
    print("Adding motifs...")
    object <- add_jaspar2024_motifs(object)

    # find differentially accessible peaks for each ident
    idents <- levels(object)
    print("Finding differentially accessible peaks...")
    for (ident in idents){
        # saved in object@misc$diff_accessible_peaks[[paste0(as.character(ident1), "_vs_", as.character(ident2))]]
        object <- diff_peaks_for_ident(object, ident1 = ident, save_in_misc = TRUE)
    }

    # find enriched motifs
    # saved in object@misc$enriched_motifs[[comparison]] with comparison in names(object@misc$diff_accessible_peaks)
    print("Finding enriched motifs...")
    object <- find_enriched_motifs(object, save_in_misc = TRUE)

    # run chromVAR
    # object[["chromvar"]]@data      # motif deviation scores (motifs × cells)
    # object[["chromvar"]]@meta.features  # motif-level metadata, if present
    print("Running chromVAR...")
    object <- RunChromVAR(object = object, genome = BSgenome.Hsapiens.UCSC.hg38, assay = "peaks")

    # find differential motif activity for each cluster
    print("Finding differential motif activity...")
    idents <- levels(object)
    for (ident in idents){
        # save in object@misc$diff_motif_activity[[paste0(as.character(ident1), "_vs_", as.character(ident2))]]
        object <- diff_motif_activity_for_ident(object, ident1 = ident, save_in_misc = TRUE)
    }

    # compute Tn5 insertion bias for footprinting, is saved in object[["peaks"]]@bias
    print("Computing Tn5 insertion bias...")
    object <- compute_Tn5_insertion_Bias(object, genome = BSgenome.Hsapiens.UCSC.hg38, assay = "peaks")

    saveRDS(object, out_path)
    invisible(NULL)
}



# Run main() only when the script is executed, not when sourced
if (sys.nframe() == 0) {
    args <- commandArgs(trailingOnly = TRUE)

    get_arg_val <- function(args, flags) {
        val <- NULL

        # handle:  --flag value   or   -f value
        for (flag in flags) {
            if (flag %in% args) {
                idx <- match(flag, args)
                if (idx < length(args)) {
                    val <- args[idx + 1]
                }
            }
        }

        # handle: --flag=value   or   -f=value
        if (is.null(val)) {
            for (flag in flags) {
                # escape dashes for regex
                pattern <- paste0("^", gsub("-", "\\\\-", flag), "=")
                eq_arg <- grep(pattern, args, value = TRUE)
                if (length(eq_arg) >= 1L) {
                    val <- sub(pattern, "", eq_arg[1L])
                    break
                }
            }
        }

        return(val)
    }

    in_path  <- get_arg_val(args, c("--in_path", "--in-path", "-i"))
    out_path <- get_arg_val(args, c("--out_path", "--out-path", "-o"))

    if (is.null(in_path) || is.null(out_path)) {
        stop(
            paste(
                "Usage:",
                "Rscript global_motif_analysis.R --in_path /path/in.rds --out_path /path/out.rds",
                "or:",
                "Rscript global_motif_analysis.R -i /path/in.rds -o /path/out.rds",
                sep = "\n"
            ),
            call. = FALSE
        )
    }

    main(in_path = in_path, out_path = out_path)
}

