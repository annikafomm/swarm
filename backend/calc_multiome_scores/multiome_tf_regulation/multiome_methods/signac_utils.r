
# Find the expected number of insertions over a set of genomic regions
# given the DNA sequences and the insertion bias of Tn5 for the experiment
# @param dna.sequence A set of DNA sequences
# @param bias A vector describing Tn5 insertion frequency at each hexamer
# @param verbose Display messages
#' @importFrom Matrix sparseMatrix
#' @importFrom Matrix colSums
#' @importFrom IRanges width narrow
FindExpectedInsertions <- function(dna.sequence, bias, verbose = TRUE) {
  if (verbose) {
    message("Computing base composition at motif sites")
  }
  total.hexamer.positions <- width(x = dna.sequence)[[1]] - 6
  hex.key <- seq_along(along.with = bias)
  names(hex.key) <- names(bias)

  # x is the hexamer frequency
  x <- vector(
    mode = "numeric",
    length = length(x = bias) * total.hexamer.positions
  )
  # i is the hexamers sequence
  i <- vector(
    mode = "numeric",
    length = length(x = bias) * total.hexamer.positions
  )
  # j is the base position
  j <- vector(
    mode = "numeric",
    length = length(x = bias) * total.hexamer.positions
  )
  current.pos <- 1

  for (jj in seq_len(length.out = total.hexamer.positions)) {
    # resize dna string set
    resized <- narrow(x = dna.sequence, start = jj, width = 6)
    resized <- as.character(x = resized)
    # need to remove any that contain N
    resized <- resized[!grepl(pattern = "N", x = resized)]
    # count
    frequencies <- table(resized)
    end.pos <- current.pos + length(x = frequencies) - 1
    # append
    x[current.pos:end.pos] <- as.numeric(x = frequencies)
    j[current.pos:end.pos] <- jj
    
    # remove frequencies not present in hex.key
    frequencies <- frequencies[names(x = frequencies) %in% names(x = hex.key)]
    
    i[current.pos:end.pos] <- as.vector(x = hex.key[names(x = frequencies)])
    # shift current position
    current.pos <- end.pos + 1
  }
  # trim vectors
  x <- x[1:(current.pos - 1)]
  i <- i[1:(current.pos - 1)]
  j <- j[1:(current.pos - 1)]

  # construct matrix
  hexamer.matrix <- sparseMatrix(
    i = i, j = j, x = x,
    dims = c(length(x = hex.key), total.hexamer.positions)
  )
  rownames(hexamer.matrix) <- names(x = hex.key)
  colnames(hexamer.matrix) <- seq_len(length.out = total.hexamer.positions)
  hexamer.matrix <- as.matrix(x = hexamer.matrix)

  if (verbose) {
    message("Computing expected Tn5 insertions per base")
  }
  # ensure correct order
  hexamer.matrix <- hexamer.matrix[names(x = bias), ]
  expected.insertions <- as.vector(
    x = crossprod(x = hexamer.matrix, y = as.matrix(x = bias))
  )

  # normalize expected by dividing by flanks
  # TODO use BackgroundMeanNorm function here
  flanks <- mean(
    x = c(expected.insertions[1:50],
          expected.insertions[
            (total.hexamer.positions - 50):total.hexamer.positions
            ])
  )
  expected.insertions <- expected.insertions / flanks
  return(expected.insertions)
}

# Create cut site pileup matrix
#
# For a set of aligned genomic ranges, find the total number of
# integration sites per cell per base.
#
# @param object A Seurat object
# @param regions A GRanges object
# @param assay Name of the assay to use
# @param cells Which cells to include. If NULL, use all cells
# @param verbose Display messages
#' @importMethodsFrom GenomicRanges strand
CreateRegionPileupMatrix <- function(
  object,
  regions,
  assay = NULL,
  cells = NULL,
  verbose = TRUE
) {
  if (length(x = regions) == 0) {
    stop("No regions supplied")
  }
  # split into strands
  on_plus <- strand(x = regions) == "+" | strand(x = regions) == "*"
  plus.strand <- regions[on_plus, ]
  minus.strand <- regions[!on_plus, ]

  # get cut matrices for each strand
  if (verbose) {
    message("Finding + strand cut sites")
  }
  cut.matrix.plus <- MultiRegionCutMatrix(
    regions = plus.strand,
    object = object,
    assay = assay,
    cells = cells,
    verbose = FALSE
  )
  if (verbose) {
    message("Finding - strand cut sites")
  }
  cut.matrix.minus <- MultiRegionCutMatrix(
    regions = minus.strand,
    object = object,
    assay = assay,
    cells = cells,
    verbose = FALSE
  )

  # reverse minus strand and add together
  if (is.null(x = cut.matrix.plus)) {
    full.matrix <- cut.matrix.minus[, rev(x = colnames(x = cut.matrix.minus))]
  } else if (is.null(x = cut.matrix.minus)) {
    full.matrix <- cut.matrix.plus
  } else {
    full.matrix <- cut.matrix.plus + cut.matrix.minus[, rev(
      x = colnames(x = cut.matrix.minus)
    )]
  }
  # rename so 0 is center
  region.width <- width(x = regions)[[1]]
  midpoint <- round(x = (region.width / 2))
  colnames(full.matrix) <- seq_len(length.out = region.width) - midpoint
  return(full.matrix)
}



# Generate cut matrix for many regions
#
# Run CutMatrix on multiple regions and add them together.
# Assumes regions are pre-aligned.
#
# @param object A Seurat object
# @param regions A set of GRanges
# @param group.by Name of grouping variable to use
# @param fragments A list of Fragment objects
# @param assay Name of the assay to use
# @param cells Vector of cells to include
# @param verbose Display messages
#' @importFrom Rsamtools TabixFile seqnamesTabix
#' @importFrom SeuratObject DefaultAssay
#' @importFrom GenomeInfoDb keepSeqlevels
MultiRegionCutMatrix <- function(
  object,
  regions,
  group.by = NULL,
  fragments = NULL,
  assay = NULL,
  cells = NULL,
  verbose = FALSE
) {
  if (inherits(x = object, what = "Seurat")) {
    assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
    object <- object[[assay]]
  }
  fragments <- SetIfNull(x = fragments, y = Fragments(object = object))
  res <- list()
  if (length(x = fragments) == 0) {
    stop("No fragment files present in assay")
  }
  for (i in seq_along(along.with = fragments)) {
    frag.path <- GetFragmentData(object = fragments[[i]], slot = "path")
    cellmap <- GetFragmentData(object = fragments[[i]], slot = "cells")
    if (is.null(x = cellmap)) {
      cellmap <- colnames(x = object)
      names(x = cellmap) <- cellmap
    }
    tabix.file <- TabixFile(
      file = frag.path,
      index = GetIndexFile(fragment = frag.path, verbose = FALSE)
    )
    open(con = tabix.file)
    # remove regions that aren't in the fragment file
    common.seqlevels <- intersect(
      x = seqlevels(x = regions),
      y = seqnamesTabix(file = tabix.file)
    )
    regions <- keepSeqlevels(
      x = regions,
      value = common.seqlevels,
      pruning.mode = "coarse"
    )
    cm <- SingleFileCutMatrix(
      cellmap = cellmap,
      tabix.file = tabix.file,
      region = regions,
      verbose = verbose
    )
    close(con = tabix.file)
    res[[i]] <- cm
  }
  # each matrix contains data for different cells at same positions
  # bind all matrices together
  res <- do.call(what = rbind, args = res)
  return(res)
}

# Set a default value if an object is null
#
# @param x An object to set if it's null
# @param y The value to provide if x is null
# @return Returns y if x is null, otherwise returns x.
SetIfNull <- function(x, y) {
  if (is.null(x = x)) {
    return(y)
  } else {
    return(x)
  }
}

# Get index file path
# checks csi and tbi extensions
# tbi used if both present
# @param fragment fragment file path
# @param verbose display messages
# @return returns the index file path
GetIndexFile <- function(fragment, verbose = TRUE) {
  is.remote <- isRemote(x = fragment)
  index.filepaths <- c(paste0(fragment, ".tbi"),
                       paste0(fragment, ".csi"))
  index.file <- index.filepaths[file.exists(index.filepaths)]
  if (length(x = index.file) == 0 & !is.remote) {
    stop("Fragment file is not indexed.")
  } else if(length(x = index.file) == 0) {
    if (verbose) {
      message("Fragment file is on a remote server")
    }
    index.file = paste0(fragment, ".tbi")
  } else if (length(x = index.file) == 2) {
    if (verbose) {
      message("TBI and CSI index both present, using TBI index")
    }
    index.file <- index.file[1]
  }
  return(index.file)
}

# Check if path is remote
# @param x path/s to check
isRemote <- function(x) {
  return(grepl(pattern = "^http|^ftp", x = x))
}

# Generate matrix of integration sites
#
# Generates a cell-by-position matrix of Tn5 integration sites
# centered on a given region (usually a DNA sequence motif). This
# matrix can be used for downstream footprinting analysis.
#
# @param cellmap A mapping of cell names in the fragment file to cell names in
# the Seurat object. Should be a named vector where each element is a cell name
# that appears in the fragment file and the name of each element is the
# name of the cell in the Seurat object.
# @param region A set of GRanges containing the regions of interest
# @param cells Which cells to include in the matrix. If NULL, use all cells in
# the cellmap
# @param tabix.file A \code{\link[Rsamtools]{TabixFile}} object.
# @param verbose Display messages
#' @importFrom Matrix sparseMatrix
#' @importFrom Rsamtools TabixFile
#' @importMethodsFrom GenomicRanges width start end
# @return Returns a sparse matrix
SingleFileCutMatrix <- function(
  cellmap,
  region,
  cells = NULL,
  tabix.file,
  verbose = TRUE
) {
  # if multiple regions supplied, must be the same width
  cells <- SetIfNull(x = cells, y = names(x = cellmap))
  if (length(x = region) == 0) {
    return(NULL)
  }
  fragments <- GetReadsInRegion(
    region = region,
    cellmap = cellmap,
    cells = cells,
    tabix.file = tabix.file,
    verbose = verbose
  )
  start.lookup <- start(x = region)
  names(start.lookup) <- seq_along(region)
  # if there are no reads in the region
  # create an empty matrix of the correct dimension
  if (nrow(x = fragments) == 0) {
    cut.matrix <- sparseMatrix(
      i = NULL,
      j = NULL,
      dims = c(length(x = cells), width(x = region)[[1]])
    )
  } else {
    fragstarts <- start.lookup[fragments$ident] + 1
    cut.df <- data.frame(
      position = c(fragments$start, fragments$end) - fragstarts,
      cell = c(fragments$cell, fragments$cell),
      stringsAsFactors = FALSE
    )
   cut.df <- cut.df[
      (cut.df$position > 0) & (cut.df$position <= width(x = region)[[1]]),
      ]
    cell.vector <- seq_along(along.with = cells)
    names(x = cell.vector) <- cells
    cell.matrix.info <- cell.vector[cut.df$cell]
    cut.matrix <- sparseMatrix(
      i = cell.matrix.info,
      j = cut.df$position,
      x = 1,
      dims = c(length(x = cells), width(x = region)[[1]])
    )
  }
  rownames(x = cut.matrix) <- cells
  colnames(x = cut.matrix) <- seq_len(width(x = region)[[1]])
  return(cut.matrix)
}

# GetReadsInRegion
#
# Extract reads for each cell within a given genomic region or set of regions
#
# @param cellmap A mapping of cell names in the fragment file to cell names in
# the Seurat object. Should be a named vector where each element is a cell name
# that appears in the fragment file and the name of each element is the
# name of the cell in the Seurat object.
# @param region A genomic region, specified as a string in the format
# 'chr:start-end'. Can be a vector of regions.
# @param tabix.file A TabixFile object.
# @param cells Cells to include. Default is all cells present in the object.
# @param verbose Display messages
# @param ... Additional arguments passed to \code{\link{StringToGRanges}}
#
#' @importFrom Rsamtools TabixFile scanTabix
#' @importFrom SeuratObject Idents
#' @importFrom fastmatch fmatch
#
# @return Returns a data frame
GetReadsInRegion <- function(
  cellmap,
  region,
  tabix.file,
  cells = NULL,
  verbose = TRUE,
  ...
) {
  file.to.object <- names(x = cellmap)
  names(x = file.to.object) <- cellmap

  if (verbose) {
    message("Extracting reads in requested region")
  }
  if (!is(object = region, class2 = "GRanges")) {
    region <- StringToGRanges(regions = region, ...)
  }
  # remove regions that aren't in the fragment file
  common.seqlevels <- intersect(
    x = seqlevels(x = region),
    y = seqnamesTabix(file = tabix.file)
  )
  if (length(common.seqlevels) != 0) {
    region <- keepSeqlevels(
      x = region,
      value = common.seqlevels,
      pruning.mode = "coarse"
    )
    reads <- scanTabix(file = tabix.file, param = region)
    reads <- TabixOutputToDataFrame(reads = reads)
    reads <- reads[
      fmatch(x = reads$cell, table = cellmap, nomatch = 0L) > 0,
    ]
    # convert cell names to match names in object
    reads$cell <- file.to.object[reads$cell]
    if (!is.null(x = cells)) {
      reads <- reads[reads$cell %in% cells, ]
    }
   if (nrow(reads) == 0) {
      reads$ident <- integer()
      reads$length <- numeric()
      return(reads)
    }
    reads$length <- reads$end - reads$start
  } else {
    reads <- data.frame(
      "chr" = character(),
      "start" = numeric(),
      "end" = numeric(),
      "cell" = character(),
      "count" = numeric(),
      "ident" = integer(),
      "length" = numeric()
    )
  }
  return(reads)
}


# TabixOutputToDataFrame
#
# Create a single dataframe from list of character vectors
#
# @param reads List of character vectors (the output of \code{\link{scanTabix}})
# @param record.ident Add a column recording which region the reads overlapped
# with
#' @importFrom stringi stri_split_fixed
#' @importFrom S4Vectors elementNROWS
TabixOutputToDataFrame <- function(reads, record.ident = TRUE) {
  if (record.ident) {
    nrep <- elementNROWS(x = reads)
  }
  original_names = names(reads)
  reads <- unlist(x = reads, use.names = FALSE)
  if (length(x = reads) == 0 | is.null(x = original_names)) {
    df <- data.frame(
      "chr" = "",
      "start" = "",
      "end" = "",
      "cell" = "",
      "count" = ""
    )
    df <- df[-1, ]
    return(df)
  }
  reads <- stri_split_fixed(str = reads, pattern = "\t")
  n <- length(x = reads[[1]])
  unlisted <- unlist(x = reads)
  e1 <- unlisted[n * (seq_along(along.with = reads)) - (n - 1)]
  e2 <- as.numeric(x = unlisted[n * (seq_along(along.with = reads)) - (n - 2)])
  e3 <- as.numeric(x = unlisted[n * (seq_along(along.with = reads)) - (n - 3)])
  e4 <- unlisted[n * (seq_along(along.with = reads)) - (n - 4)]
  e5 <- as.numeric(x = unlisted[n * (seq_along(along.with = reads)) - (n - 5)])
  df <- data.frame(
    "chr" = e1,
    "start" = e2,
    "end" = e3,
    "cell" = e4,
    "count" = e5,
    stringsAsFactors = FALSE,
    check.rows = FALSE,
    check.names = FALSE
  )
  if (record.ident) {
    df$ident <- rep(x = seq_along(along.with = nrep), nrep)
  }
  return(df)
}