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
    spot_mean <- Diagonal(x = 1 / pmax(w_spot, 1e-12)) %*% spot_sum
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

