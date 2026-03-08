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
    out_path,
    ATAC_assay = "peaks",
    species = "Homo sapiens",
    genome = BSgenome.Hsapiens.UCSC.hg38
    ){
    if (is.null(object[[ATAC_assay]]@motifs)){
        jaspar <- JASPAR2024()
        sq24 <- RSQLite::dbConnect(RSQLite::SQLite(), db(jaspar))
        pfm <- TFBSTools::getMatrixSet(sq24, list(species = species, collection = "CORE", all_versions = FALSE, matrixtype = "PFM"))
        object <- AddMotifs(object, genome = genome, pfm = pfm, assay=ATAC_assay)
        motif_to_tf <- data.frame(
          motif_id = vapply(pfm, TFBSTools::ID, character(1)),
          TF       = vapply(pfm, TFBSTools::name, character(1)),
          stringsAsFactors = FALSE
        )
        # save motif to tf data.frame
        write.csv(motif_to_tf, file = file.path(out_path, "motif_to_tf.csv"), row.names = FALSE)
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


#' Project chromVAR motif deviation scores from dissociated cells onto spatial spots
#'
#' ChromVAR deviation scores computed per dissociated cell are projected onto
#' spatial spots via a weighted mean using the Tangram mapping matrix.
#' The resulting spot-level scores are stored as a new \code{chromvar} assay
#' inside \code{spot_obj} and active idents are set to \code{clustering_var}.
#'
#' @param object_dissociated A Seurat object with a chromVAR assay
#'   (motifs \eqn{\times} cells deviation scores).
#' @param spot_obj A spot-level Seurat object created with
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
    #spot_obj,
    M,
    chromvar_assay = "chromvar",
    clustering_var = "leiden"
) {
    # MERIT 
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

    # Add metadata (now aligned)
    spot_obj <- AddMetaData(spot_obj, spot_meta)
    return(spot_obj)

    # # Both spot_obj and the Tangram M columns may use different barcode
    # # suffixes (".1" vs "-1" is common in read.csv output vs 10x format).
    # # Normalise all barcodes to the "-N" form before any matching.
    # fix_barcode <- function(x) sub("\\.(\\d+)$", "-\\1", x)

    # # motifs × cells deviation matrix
    # cv_mat      <- GetAssayData(object_dissociated, assay = chromvar_assay, layer = "data")

    # # align cells between cv_mat and M
    # shared_cells <- intersect(colnames(cv_mat), rownames(M))
    # if (length(shared_cells) == 0L)
    #     stop("project_chromvar_to_spots: no shared cell barcodes between chromVAR assay and Tangram M matrix.")
    # cv_mat <- cv_mat[, shared_cells, drop = FALSE]
    # M_mat  <- as.matrix(M[shared_cells, , drop = FALSE])

    # # Normalise M column names (spots)
    # colnames(M_mat) <- fix_barcode(colnames(M_mat))

    # # weighted mean: (motifs × cells) %*% (cells × spots) then normalise
    # M_sp      <- Matrix::Matrix(M_mat, sparse = TRUE)
    # spot_sum  <- cv_mat %*% M_sp                               # motifs × spots
    # w_spot    <- Matrix::colSums(M_sp)
    # spot_mean <- spot_sum %*% Matrix::Diagonal(x = 1 / pmax(w_spot, 1e-12))
    # colnames(spot_mean) <- colnames(M_sp)                      # already fixed

    # # Normalise spot_obj barcodes for matching (Seurat may ignore colnames<-)
    # spot_obj_cells_raw   <- colnames(spot_obj)                 # original names
    # spot_obj_cells_fixed <- fix_barcode(spot_obj_cells_raw)

    # common_fixed <- intersect(colnames(spot_mean), spot_obj_cells_fixed)
    # if (length(common_fixed) == 0L) {
    #     message("spot_mean columns (first 5): ", paste(head(colnames(spot_mean), 5), collapse = ", "))
    #     message("spot_obj cells (first 5):    ", paste(head(spot_obj_cells_fixed, 5), collapse = ", "))
    #     stop("project_chromvar_to_spots: no overlapping spot barcodes between projected scores and spot_obj.")
    # }

    # # Subset spot_mean to common spots
    # spot_mean <- spot_mean[, common_fixed, drop = FALSE]

    # # Subset spot_obj using the ORIGINAL (un-fixed) names that Seurat knows
    # orig_idx   <- match(common_fixed, spot_obj_cells_fixed)
    # spot_obj   <- spot_obj[, spot_obj_cells_raw[orig_idx]]

    # # Rename spot_mean columns to match whatever Seurat stored
    # colnames(spot_mean) <- colnames(spot_obj)

    # # store projected scores as a chromvar assay in spot_obj
    # spot_obj[[chromvar_assay]] <- CreateAssayObject(data = spot_mean)

    # Seurat::Idents(spot_obj) <- clustering_var
    # return(spot_obj)
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

#' Project dissociated-cell ATAC footprint profiles onto spatial spots
#'
#' This function projects transcription factor footprint position-enrichment
#' profiles computed on dissociated scATAC / multiome data (via Signac's
#' \code{Footprint()}) onto spatial transcriptomics spots using a Tangram
#' cell-to-spot mapping matrix.
#'
#' The projection is performed as a weighted mean of per-cell footprint
#' profiles, where the weights correspond to Tangram mapping probabilities.
#' The resulting matrix has spatial spots as rows and footprint positions
#' (e.g. -250 to +250 bp around motif centers) as columns, with the original
#' \code{expected} and \code{motif} rows appended for compatibility with
#' Signac's footprint plotting utilities.
#'
#' The returned footprint profiles represent *mixtures of
#' dissociated-cell ATAC signals*. They are not footprints computed from
#' true spatial ATAC fragments and should be interpreted accordingly.
#'
#' @param object_dissociated A \code{Seurat} object containing dissociated
#'   scATAC-seq or multiome data with a \code{ChromatinAssay}. The assay must
#'   already contain footprint position enrichment data computed using
#'   \code{Signac::Footprint()}.
#'
#' @param footprint_key Character scalar. Name of the footprint entry stored
#'   in \code{object_dissociated[[assay]]@positionEnrichment}, e.g.
#'   \code{"footprint_MA0084.2"}.
#'
#' @param M A numeric matrix or sparse matrix of dimensions
#'   \code{(n_cells × n_spots)} containing Tangram mapping probabilities.
#'   Row names must correspond to dissociated cell barcodes (matching
#'   \code{colnames(object_dissociated)}), and column names must correspond
#'   to spatial spot identifiers.
#'
#' @param assay Character scalar. Name of the ChromatinAssay in
#'   \code{object_dissociated} from which footprint data should be read.
#'   Defaults to \code{"peaks"}.
#'
#' @return A \code{dgCMatrix} with dimensions
#'   \code{(n_spots + 2) × n_positions}. Rows correspond to spatial spot IDs,
#'   followed by the \code{"expected"} and \code{"motif"} rows copied from the
#'   original footprint matrix. Columns correspond to positions relative to
#'   the motif center.
#'
#' @details
#' The projection is computed as:
#' \deqn{
#'   F_{spot} = D^{-1} (M^\top F_{cell})
#' }
#' where \eqn{F_{cell}} is the dissociated-cell footprint matrix,
#' \eqn{M} is the Tangram mapping matrix, and \eqn{D} is a diagonal matrix
#' of per-spot weight sums used for normalization.
#'
#' @seealso \code{\link[Signac]{Footprint}}, \code{\link[Signac]{PlotFootprint}}
#'
#' @examples
#' \dontrun{
#' fp_spatial <- footprints_dissociated2spatial(
#'   object_dissociated = processed_data,
#'   footprint_key = "footprint_7_vs_rest_MA0084.2",
#'   M = tangram_map
#' )
#' }
#' Construct a spatial Seurat object by projecting dissociated ATAC peak counts
#'
#' This function creates a new \code{Seurat} object whose "cells" correspond to
#' spatial transcriptomics spots, by projecting dissociated scATAC-seq or multiome
#' peak counts onto spatial spots using a Tangram cell-to-spot mapping matrix.
#'
#' Peak counts are projected as a weighted sum of dissociated-cell peak counts,
#' where weights are given by Tangram mapping probabilities. The resulting object
#' contains a \code{ChromatinAssay} with real genomic peak ranges but *projected*
#' (non-integer) counts.
#'
#' The returned object does *not* contain spatial ATAC fragments.
#' As a consequence, operations requiring fragment-level information (e.g.
#' \code{Signac::Footprint()}, Tn5 insertion bias computation) are not valid on
#' this object. The object is intended for visualization, aggregation, and
#' downstream integration with precomputed footprint or motif data.
#'
#' @param object_dissociated A \code{Seurat} object containing dissociated
#'   scATAC-seq or multiome data with a \code{ChromatinAssay}. The assay must
#'   contain peak-level counts and genomic ranges.
#'
#' @param M A numeric matrix or sparse matrix of dimensions
#'   \code{(n_cells × n_spots)} containing Tangram mapping probabilities.
#'   Row names must correspond to dissociated cell barcodes (matching
#'   \code{colnames(object_dissociated)}), and column names must correspond
#'   to spatial spot identifiers.
#'
#' @param spot_meta A \code{data.frame} containing spatial spot metadata
#'   (e.g. spatial clusters, tissue regions). Row names must exactly match
#'   the spatial spot identifiers (column names of \code{M}).
#'
#' @param assay Character scalar. Name of the ChromatinAssay in
#'   \code{object_dissociated} to project. Defaults to \code{"peaks"}.
#'
#' @param slot Character scalar. Slot of the ChromatinAssay to use as the
#'   source of peak-level values. Defaults to \code{"counts"}.
#'
#' @return A \code{Seurat} object whose cells correspond to spatial spots and
#'   whose default assay is a \code{ChromatinAssay} containing projected
#'   peak counts and the original genomic peak ranges.
#'
#' @details
#' The projection is computed as:
#' \deqn{
#'   C_{spot} = C_{cell} \times M
#' }
#' where \eqn{C_{cell}} is the peak-by-cell matrix from the dissociated object
#' and \eqn{M} is the Tangram cell-to-spot mapping matrix.
#'
#' The resulting counts are continuous-valued and should not be interpreted
#' as observed fragment counts.
#'
#' @seealso \code{\link[Signac]{CreateChromatinAssay}},
#'   \code{\link[Seurat]{CreateSeuratObject}},
#'   \code{\link[Signac]{Footprint}}
#'
#' @examples
#' \dontrun{
#' spot_obj <- seuratObj_dissociated2spatial(
#'   object_dissociated = processed_data,
#'   M = tangram_map,
#'   spot_meta = spatial_metadata
#' )
#' }
seuratObj_dissociated2spatial <- function(
    object_dissociated,
    M, # tangram map eg M = read.csv("../swarm/backend/uploads/job_1768478211396_76193c44-3781-4b8b-b9ef-fc7702b4a617/plasmidpoop/tangram_map.csv",row.names=1), adata_map.X
    spot_meta, # <-read.csv("../swarm/backend/uploads/job_1768478211396_76193c44-3781-4b8b-b9ef-fc7702b4a617/plasmidpoop/meta_spatial.csv"), adata_map.var
    assay = "peaks",
    slot = "counts"
){
    C <- GetAssayData(object_dissociated, assay = assay, slot = slot)
    # force to 2D numeric (and to sparse for speed)
    M_sp <- Matrix::Matrix(as.matrix(M), sparse = TRUE)

    C_spot <- C %*% M_sp

    ranges <- granges(object_dissociated[[assay]])  # peak GRanges from your multiome object

    chrom_spot <- CreateChromatinAssay(
        counts = C_spot,
        ranges = ranges,
        genome = genome(object_dissociated[[assay]]) # or "hg38"
    )

    spot_obj <- CreateSeuratObject(counts = chrom_spot, assay = assay)
    DefaultAssay(spot_obj) <- assay

    spot_obj <- AddMetaData(spot_obj, spot_meta)
    fix_barcode <- function(x) sub("\\.(\\d+)$", "-\\1", x)
    # 1) fix spot IDs consistently
    colnames(spot_obj) <- fix_barcode(colnames(spot_obj))
    return(spot_obj)

}

to_dgC <- function(x) {
  if (inherits(x, "dgCMatrix")) return(x)
  # works for data.frame / matrix / dense Matrix
  Matrix::Matrix(as.matrix(x), sparse = TRUE)
}

footprints_dissociated2spatial <- function(
    object_dissociated,
    footprint_key,
    M, # tangram map eg M = read.csv("../swarm/backend/uploads/job_1768478211396_76193c44-3781-4b8b-b9ef-fc7702b4a617/plasmidpoop/tangram_map.csv",row.names=1)
    assay = "peaks"
    ){

    # 1) get the footprint matrix (cells + expected + motif) x positions
    fp_all <- object_dissociated[[assay]]@positionEnrichment[[footprint_key]]

    # 2) remove motif data and bias
    cells <- colnames(object_dissociated)
    fp <- fp_all[intersect(rownames(fp_all), cells), , drop = FALSE]   # cells x positions
    # Make sure ordering matches
    fp <- fp[rownames(M), , drop = FALSE]

    # 4) weighted SUM to spots: (spots x cells) %*% (cells x positions) = spots x positions
    M = Matrix::Matrix(as.matrix(M), sparse = TRUE)
    spot_sum <- t(M) %*% fp

    # normalize
    w_spot <- Matrix::rowSums(t(M))  # length = n_spots
    spot_mean <- spot_sum
    spot_mean <- Matrix::Diagonal(x = 1 / pmax(w_spot, 1e-12)) %*% spot_sum
    spot_mean_sp <- spot_mean
    if (!inherits(spot_mean_sp, "dgCMatrix")) {
        spot_mean_sp <- Matrix::Matrix(as.matrix(spot_mean_sp), sparse = TRUE)
    }
    rownames(spot_mean_sp) <- rownames(spot_sum)

    # 5 add bias and motif back
    expected <- fp_all["expected", , drop = FALSE]
    motif <- fp_all["motif", , drop = FALSE]
    fp_spots_all <- rbind(spot_mean_sp, expected, motif)
    return(fp_spots_all)
    }


#' Compute and plot spatially projected ATAC footprints for a single motif, and
#' save updated Seurat objects to disk
#'
#' This function loads a dissociated (multiome/scATAC) Seurat object and a
#' spot-level Seurat object from disk, computes a motif footprint in the
#' dissociated object via \code{Signac::Footprint()}, projects the footprint to
#' spatial spots using a Tangram cell-to-spot mapping matrix, stores the
#' projected footprint matrix in the spot object, and returns a
#' \code{Signac::PlotFootprint()} plot.
#'
#' Unlike the in-memory workflow, this path-based variant **persists changes**
#' by writing the updated Seurat objects back to disk.
#'
#' The spatial footprints shown are **projected mixtures** of dissociated-cell
#' ATAC footprints (weighted by Tangram probabilities). They do not represent
#' footprints computed from true spatial ATAC fragments.
#'
#' @param object_dissociated_path Character scalar. Path to an RDS file
#'   containing a \code{Seurat} object with dissociated scATAC-seq or multiome
#'   data (must be compatible with \code{Signac::Footprint()}).
#'
#' @param M_path Character scalar. Path to a CSV file containing the Tangram
#'   mapping matrix (\code{n_cells × n_spots}). Row names must match dissociated
#'   cell barcodes; column names must match spatial spot IDs.
#'
#' @param spot_obj_path Character scalar. Path to an RDS file containing a
#'   spot-level \code{Seurat} object (cells = spatial spots). Must contain an
#'   assay named \code{assay} and a metadata column \code{clustering_var}.
#'
#' @param motif_name Character scalar. Motif identifier/name passed to
#'   \code{Signac::Footprint(motif.name = ...)} (e.g. \code{"MA0645.2"}).
#'
#' @param assay Character scalar. Name of the \code{ChromatinAssay} used for
#'   footprint computation and for storing projected footprints in the spot
#'   object. Defaults to \code{"peaks"}.
#'
#' @param clustering_var Character scalar. Column name in
#'   \code{spot_obj@meta.data} used to group spots in the plot. Defaults to
#'   \code{"leiden"}.
#'
#' @param object_dissociated_out_path Character scalar. Output path for saving
#'   the updated dissociated Seurat object. Defaults to
#'   \code{object_dissociated_path} (overwrite).
#'
#' @param spot_obj_out_path Character scalar. Output path for saving the updated
#'   spot-level Seurat object. Defaults to \code{spot_obj_path} (overwrite).
#'
#' @param overwrite Logical. Whether to overwrite existing output files.
#'   Defaults to \code{TRUE}.
#'
#' @return A \code{ggplot} object returned by \code{Signac::PlotFootprint()}.
#'
#' @seealso \code{\link[Signac]{Footprint}}, \code{\link[Signac]{PlotFootprint}},
#'   \code{footprints_dissociated2spatial}
plot_footprint_for_motif <- function(
  M_path, # eg adata_map.X.csv
  spot_obj_path,
  motif_name,
  object_dissociated = NULL,
  object_dissociated_path = NULL,
  assay = "peaks",
  clustering_var = "leiden",
  object_dissociated_out_path = object_dissociated_path,
  spot_obj_out_path = spot_obj_path,
  plot_out_path = NULL,
  overwrite = TRUE
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
  M <- read.csv(M_path, row.names = 1, check.names = FALSE)
  spot_obj <- readRDS(spot_obj_path)

  key <- sprintf("footprint_%s", paste(motif_name))

  # Compute footprint in dissociated object
  object_dissociated <- Signac::Footprint(
    object = object_dissociated,
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

  # Store in spot object
  spot_obj[[assay]]@positionEnrichment[[key]] <- fp_spots_all
  Seurat::Idents(spot_obj) <- clustering_var

  # Persist updated objects to disk
  saveRDS(object_dissociated, file = object_dissociated_out_path)
  saveRDS(spot_obj, file = spot_obj_out_path)

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

  return(p)
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
    object <- add_jaspar2024_motifs(object, out_path)

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




# =============================================================================
# Patch: Spatial projected chromVAR + footprinting helpers (ported from
#        test_diff_motif_activityR_healthy*.ipynb)
#        (Placed at end so these definitions override earlier ones.)
# =============================================================================

# Global helper: harmonize Visium spot barcodes that often appear as "AAAC.1"
# in CSV exports to the canonical "AAAC-1"
fix_barcode <- function(x) {
  sub("\\.(\\d+)$", "-\\1", x)
}

# Seurat v4/v5 compatibility helpers
get_assay_data_compat <- function(object, assay, layer = "data", slot_fallback = "data") {
  out <- tryCatch({
    GetAssayData(object = object, assay = assay, layer = layer)
  }, error = function(e) NULL)

  if (is.null(out)) {
    out <- GetAssayData(object = object, assay = assay, slot = slot_fallback)
  }
  return(out)
}

set_assay_data_compat <- function(object, assay, layer = "data", new.data) {
  # Seurat v5 prefers layer; v4 uses slot.
  object <- tryCatch({
    SetAssayData(object = object, assay = assay, layer = layer, new.data = new.data)
  }, error = function(e) {
    SetAssayData(object = object, assay = assay, slot = layer, new.data = new.data)
  })
  return(object)
}

standardize_spot_meta <- function(spot_meta,
                                  spot_id_col = "X",
                                  fix_spot_ids = TRUE) {
  if (is.null(spot_meta)) return(NULL)

  spot_meta <- as.data.frame(spot_meta)

  # If rownames are missing/boring, use a column (default 'X' like read.csv output)
  if (is.null(rownames(spot_meta)) || all(rownames(spot_meta) %in% c("", seq_len(nrow(spot_meta))))) {
    if (!is.null(spot_id_col) && spot_id_col %in% colnames(spot_meta)) {
      rownames(spot_meta) <- as.character(spot_meta[[spot_id_col]])
    }
  }

  if (fix_spot_ids) {
    rownames(spot_meta) <- fix_barcode(rownames(spot_meta))
  }

  return(spot_meta)
}

align_M_to_cells <- function(M, cell_ids) {
  # M is cells x spots
  if (is.data.frame(M)) M <- as.matrix(M)

  if (is.null(rownames(M))) {
    stop("Tangram map M must have rownames matching single-cell barcodes / cell IDs.")
  }

  common <- intersect(cell_ids, rownames(M))
  if (length(common) < 10) {
    stop(sprintf("Too few overlapping cells between chromVAR matrix (%d) and M rownames (%d). Overlap=%d",
                 length(cell_ids), nrow(M), length(common)))
  }

  # Reorder to match chromVAR columns
  M <- M[common, , drop = FALSE]
  return(list(M = M, common_cells = common))
}

project_matrix_to_spots <- function(feature_by_cell, M) {
  # feature_by_cell: features x cells
  # M: cells x spots
  if (is.data.frame(M)) M <- as.matrix(M)
  M_sp <- Matrix::Matrix(M, sparse = TRUE)
  feature_by_spot <- feature_by_cell %*% M_sp
  return(feature_by_spot)
}

# # ---- Differential motif activity on spatial (Visium) spots ----
# # Build a Seurat object whose "cells" are Visium spots and whose assay is projected chromVAR deviations.
# project_chromvar_to_spots <- function(object_dissociated,
#                                       M,
#                                       spot_meta = NULL,
#                                       chromvar_assay = "chromvar",
#                                       spot_id_col = "X",
#                                       fix_spot_ids = TRUE) {
#   C <- get_assay_data_compat(object_dissociated, assay = chromvar_assay, layer = "data", slot_fallback = "data")
#   # C: motifs x cells
#   aln <- align_M_to_cells(M, cell_ids = colnames(C))
#   M_use <- aln$M
#   C_use <- C[, aln$common_cells, drop = FALSE]

#   C_spot <- project_matrix_to_spots(C_use, M_use)  # motifs x spots

#   if (fix_spot_ids) {
#     colnames(C_spot) <- fix_barcode(colnames(C_spot))
#   }

#   # Ensure sparse dgCMatrix for Seurat
#   C_spot <- to_dgC(C_spot)

#   spot_obj <- CreateSeuratObject(counts = C_spot, assay = chromvar_assay)
#   spot_obj <- set_assay_data_compat(spot_obj, assay = chromvar_assay, layer = "data", new.data = C_spot)

#   # Attach metadata
#   spot_meta <- standardize_spot_meta(spot_meta, spot_id_col = spot_id_col, fix_spot_ids = fix_spot_ids)
#   if (!is.null(spot_meta)) {
#     common_spots <- intersect(colnames(spot_obj), rownames(spot_meta))
#     spot_obj <- subset(spot_obj, cells = common_spots)
#     spot_meta <- spot_meta[common_spots, , drop = FALSE]
#     spot_obj <- AddMetaData(spot_obj, metadata = spot_meta)
#   }

#   return(spot_obj)
# }

# Override: add logfc.threshold + only_pos default like in notebook
diff_motif_activity_for_ident <- function(object,
                                         ident1,
                                         ident2 = NULL,
                                         save_in_misc = TRUE,
                                         chromvar_assay = "chromvar",
                                         only_pos = FALSE,
                                         logfc.threshold = 0.1) {
  da <- FindMarkers(
    object,
    ident.1 = ident1,
    ident.2 = ident2,
    assay = chromvar_assay,
    slot = "data",
    only.pos = only_pos,
    logfc.threshold = logfc.threshold
  )
  if (save_in_misc) {
    if (is.null(object@misc$diff_motif_activity)) object@misc$diff_motif_activity <- list()
    comparison <- if (is.null(ident2)) {
      paste0(as.character(ident1), "_vs_all")
    } else {
      paste0(as.character(ident1), "_vs_", as.character(ident2))
    }
    object@misc$diff_motif_activity[[comparison]] <- da
  }
  return(object)
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
      only_pos = only_pos,
      logfc.threshold = logfc.threshold
    )
  }
  return(spot_obj)
}

top_motifs_per_comparison <- function(spot_obj, n = 3) {
  if (is.null(spot_obj@misc$diff_motif_activity)) {
    return(data.frame())
  }
  out <- lapply(names(spot_obj@misc$diff_motif_activity), function(comp) {
    df <- spot_obj@misc$diff_motif_activity[[comp]]
    if (is.null(df) || nrow(df) == 0) return(NULL)
    df$motif <- rownames(df)
    df$comparison <- comp
    # Seurat FindMarkers returns 'p_val_adj' and 'avg_log2FC' or 'avg_diff' depending on version
    # Prefer avg_log2FC if present, else avg_diff.
    eff <- if ("avg_log2FC" %in% colnames(df)) "avg_log2FC" else if ("avg_diff" %in% colnames(df)) "avg_diff" else NULL
    if (!is.null(eff)) df <- df[order(df$p_val_adj, -df[[eff]]), , drop = FALSE]
    head(df[, intersect(c("comparison", "motif", "p_val_adj", "p_val", "avg_log2FC", "avg_diff"), colnames(df))], n)
  })
  out <- do.call(rbind, out)
  if (is.null(out)) out <- data.frame()
  rownames(out) <- NULL
  return(out)
}

# ---- Footprinting on projected (spot-level) chromatin assay ----
# Add a projected chromatin assay (e.g. peaks counts) to an existing spot Seurat object.
add_projected_atac_assay_to_spot_obj <- function(object_dissociated,
                                                 M,
                                                 spot_obj,
                                                 assay = "peaks",
                                                 slot = "counts",
                                                 fix_spot_ids = TRUE) {
  counts <- GetAssayData(object = object_dissociated, assay = assay, slot = slot)  # peaks x cells
  aln <- align_M_to_cells(M, cell_ids = colnames(counts))
  M_use <- aln$M
  counts_use <- counts[, aln$common_cells, drop = FALSE]

  counts_spot <- project_matrix_to_spots(counts_use, M_use)  # peaks x spots

  if (fix_spot_ids) {
    colnames(counts_spot) <- fix_barcode(colnames(counts_spot))
  }
  counts_spot <- to_dgC(counts_spot)

  # Keep only spots present in spot_obj
  common_spots <- intersect(colnames(spot_obj), colnames(counts_spot))
  spot_obj <- subset(spot_obj, cells = common_spots)
  counts_spot <- counts_spot[, common_spots, drop = FALSE]

  chrom_assay <- CreateChromatinAssay(
    counts = counts_spot,
    ranges = granges(object_dissociated[[assay]]),
    genome = genome(object_dissociated[[assay]]),
    motifs = Motifs(object_dissociated[[assay]]),
    annotation = Annotation(object_dissociated[[assay]])
  )
  spot_obj[[assay]] <- chrom_assay
  return(spot_obj)
}

# Compute, project, and plot a footprint for a single motif, working entirely in-memory.
# Saves a PDF if plot_out_path is provided.
compute_project_plot_footprint <- function(motif_name,
                                           object_dissociated,
                                           spot_obj,
                                           M,
                                           assay = "peaks",
                                           clustering_var = "cell_type",
                                           plot_out_path = NULL,
                                           overwrite = TRUE) {
  if (!clustering_var %in% colnames(spot_obj@meta.data)) {
    stop(sprintf("clustering_var='%s' not found in spot_obj@meta.data.", clustering_var))
  }
  Idents(spot_obj) <- spot_obj[[clustering_var, drop = TRUE]]

  # Compute footprint on dissociated object (stores in object_dissociated[[assay]]@positionEnrichment)
  key <- paste0("footprint_", motif_name)
  object_dissociated <- Footprint(
    object = object_dissociated,
    motif.name = motif_name,
    assay = assay,
    key = key
  )

  # Project to spots
  fp_spots_all <- footprints_dissociated2spatial(object_dissociated, footprint_key = key, M = M, assay = assay)

  # Fix spot barcode style inside the footprint matrix (only for spot rows)
  spot_rows <- setdiff(rownames(fp_spots_all), c("expected", "motif"))
  rn <- rownames(fp_spots_all)
  rn[rn %in% spot_rows] <- fix_barcode(rn[rn %in% spot_rows])
  rownames(fp_spots_all) <- rn

  # Store in spot_obj assay
  spot_obj[[assay]]@positionEnrichment[[key]] <- fp_spots_all

  p <- PlotFootprint(
    object = spot_obj,
    features = motif_name,
    assay = assay,
    group.by = clustering_var,
    key = key,
    show.expected = TRUE
  )

  if (!is.null(plot_out_path)) {
    if (file.exists(plot_out_path) && !overwrite) {
      message("plot_out_path exists and overwrite=FALSE: ", plot_out_path)
    } else {
      ggsave(filename = plot_out_path, plot = p, width = 9, height = 5)
    }
  }

  return(list(plot = p, spot_obj = spot_obj, object_dissociated = object_dissociated, key = key))
}


# Override for robustness: intersect cells between footprint matrix and M before projecting
footprints_dissociated2spatial <- function(object_dissociated,
                                          footprint_key,
                                          M,
                                          assay = "peaks") {
  if (is.data.frame(M)) M <- as.matrix(M)
  if (is.null(rownames(M))) {
    stop("Tangram map M must have rownames matching dissociated cell IDs.")
  }

  fp_all <- object_dissociated[[assay]]@positionEnrichment[[footprint_key]]
  if (is.null(fp_all)) {
    stop("No footprint found in object_dissociated[[", assay, "]]@positionEnrichment[['", footprint_key, "']].")
  }

  # Keep only real cells (exclude expected/motif rows)
  cell_ids <- colnames(object_dissociated)
  fp_cells <- fp_all[intersect(rownames(fp_all), cell_ids), , drop = FALSE]

  common <- intersect(rownames(M), rownames(fp_cells))
  if (length(common) < 10) {
    stop("Too few overlapping cells between M rownames and footprint rows. Overlap=", length(common))
  }

  fp_cells <- fp_cells[common, , drop = FALSE]
  M_use <- M[common, , drop = FALSE]
  M_use <- Matrix::Matrix(M_use, sparse = TRUE)

  spot_sum <- t(M_use) %*% fp_cells
  w_spot <- Matrix::rowSums(t(M_use))
  w_spot[w_spot == 0] <- 1
  spot_mean <- Matrix::Diagonal(x = 1 / w_spot) %*% spot_sum
  spot_mean_sp <- to_dgC(spot_mean)
  rownames(spot_mean_sp) <- rownames(spot_sum)

  expected <- fp_all["expected", , drop = FALSE]
  motif <- fp_all["motif", , drop = FALSE]
  fp_spots_all <- rbind(spot_mean_sp, expected, motif)

  return(fp_spots_all)
}
