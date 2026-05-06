# ============================================================
# Parallel-safe motif enrichment rewrite
# ============================================================

suppressPackageStartupMessages({
  library(future)
  library(future.apply)
  library(GenomicRanges)
  library(IRanges)
})

# ------------------------------------------------------------
# small helpers
# ------------------------------------------------------------

.make_gene_cluster_key <- function(gene, cluster) {
  paste(gene, cluster)
}

.peaks_to_granges <- function(peaks) {
  if (length(peaks) == 0) {
    return(GRanges())
  }
  peaks_data <- strsplit(peaks, "-")
  chromosome <- vapply(peaks_data, `[`, character(1), 1)
  start_pos  <- as.numeric(vapply(peaks_data, `[`, character(1), 2))
  end_pos    <- as.numeric(vapply(peaks_data, `[`, character(1), 3))
  GRanges(seqnames = chromosome, ranges = IRanges(start = start_pos, end = end_pos))
}

.count_motifs_in_peaks <- function(motif_data, peaks, col_names = colnames(motif_data)) {
  out <- setNames(numeric(length(col_names)), col_names)

  if (length(peaks) == 0) {
    return(out)
  }

  mat <- motif_data[peaks, , drop = FALSE]
  if (inherits(mat, "sparseMatrix")) {
    mat <- as.matrix(mat)
  }

  if (nrow(mat) > 1) {
    vals <- colSums(mat == TRUE)
  } else {
    vals <- as.numeric(mat[1, ] == TRUE)
    names(vals) <- col_names
  }

  out[names(vals)] <- as.numeric(vals)
  out
}

.safe_t_test_less <- function(x, mu) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]

  if (length(x) < 2) {
    return(c(statistic = NA_real_, p.value = 1))
  }

  res <- tryCatch(
    t.test(x = x, mu = mu, alternative = "less"),
    error = function(e) NULL
  )

  if (is.null(res)) {
    return(c(statistic = NA_real_, p.value = 1))
  }

  c(statistic = unname(res$statistic), p.value = unname(res$p.value))
}

.rbind_fill_df <- function(df1, df2) {
  if (is.null(df1)) return(df2)
  if (is.null(df2)) return(df1)

  all_cols <- union(colnames(df1), colnames(df2))
  miss1 <- setdiff(all_cols, colnames(df1))
  miss2 <- setdiff(all_cols, colnames(df2))

  for (cc in miss1) df1[[cc]] <- NA
  for (cc in miss2) df2[[cc]] <- NA

  df1 <- df1[, all_cols, drop = FALSE]
  df2 <- df2[, all_cols, drop = FALSE]

  rbind(df1, df2)
}

# top-level worker, outside calculate_enrichments()
.motif_worker <- function(
    i,
    helper.df,
    gene_peak_links,
    motif_data,
    peak_meta,
    motifs_in_proximal,
    motifs_in_distal,
    motif_granges_cache,
    bg_size
) {
  gene <- helper.df$gene[[i]]
  cluster <- helper.df$cluster[[i]]

  motif_enrichment_per_gene_df(
    gene = gene,
    cluster = cluster,
    gene_peak_links = gene_peak_links,
    motif_data = motif_data,
    peak_meta = peak_meta,
    motifs_in_proximal = motifs_in_proximal,
    motifs_in_distal = motifs_in_distal,
    motif_granges_cache = motif_granges_cache,
    bg_size = bg_size
  )
}

# ------------------------------------------------------------
# fixed version of compute_motifs_to_gene
# ------------------------------------------------------------

compute_motifs_to_gene <- function(
    object, gene, cluster,
    assay = "peaks", regulatorType = "proximal",
    gene_peak_links = NULL, overwrite = FALSE,
    test.in_promoter = FALSE, filtered = TRUE
) {
  # Computes which motifs are located in peaks associated to gene.
  # Stores GRangesList-like results in:
  #   object@misc$motifs.inProximalPeaks[[paste(gene, cluster)]]
  #   object@misc$motifs.inDistalPeaks[[paste(gene, cluster)]]

  key <- .make_gene_cluster_key(gene, cluster)

  if (is.null(object@misc$motifs.inProximalPeaks)) object@misc$motifs.inProximalPeaks <- list()
  if (is.null(object@misc$motifs.inDistalPeaks))   object@misc$motifs.inDistalPeaks   <- list()

  if (regulatorType == "proximal" &&
      !is.null(object@misc$motifs.inProximalPeaks[[key]]) &&
      !overwrite) {
    return(object)
  }

  if (regulatorType == "distal" &&
      !is.null(object@misc$motifs.inDistalPeaks[[key]]) &&
      !overwrite) {
    return(object)
  }

  if (is.null(gene_peak_links)) {
    gene_peak_links <- object@misc$peak_stats
  }

  if (test.in_promoter && regulatorType == "distal") {
    peaks_linked_to_gene <- gene_peak_links$peak[
      gene_peak_links$gene == gene &
      gene_peak_links$cluster == cluster
    ]
  } else {
    peaks_linked_to_gene <- gene_peak_links$peak[
      gene_peak_links$gene == gene &
      gene_peak_links$cluster == cluster &
      gene_peak_links$regulatorType == regulatorType
    ]
  }

  # store empty list to avoid recomputing forever
  if (length(peaks_linked_to_gene) == 0) {
    if (regulatorType == "proximal") {
      object@misc$motifs.inProximalPeaks[[key]] <- list()
    } else {
      object@misc$motifs.inDistalPeaks[[key]] <- list()
    }
    return(object)
  }

  peaks_GRanges <- .peaks_to_granges(peaks_linked_to_gene)

  motif_data <- object[[assay]]@motifs@data
  motif_count <- .count_motifs_in_peaks(motif_data, peaks_linked_to_gene)
  motifs_for_gene <- names(motif_count)[motif_count > 0]

  if (length(motifs_for_gene) == 0) {
    result_list <- list()
  } else {
    motif_GRangesList <- get_GRanges_for_motifs(object, motifs_for_gene)
    result_list <- lapply(motif_GRangesList, function(gr) {
      subsetByOverlaps(gr, peaks_GRanges)
    })
    result_list <- result_list[vapply(result_list, length, integer(1)) > 0]
  }

  if (test.in_promoter && length(result_list) > 0) {
    geneCoords <- genes(EnsDb.Hsapiens.v86, filter = GeneNameFilter(gene))
    geneStrand <- as.character(strand(geneCoords))
    geneChrom  <- paste0("chr", as.character(seqnames(geneCoords)))

    if (length(geneStrand) > 1) {
      geneCoords <- geneCoords[1]
      geneStrand <- geneStrand[1]
      geneChrom  <- geneChrom[1]
    }

    if (geneStrand == "-") {
      promStart <- end(geneCoords) + 20
      promEnd   <- end(geneCoords) + 2000
    } else {
      promEnd   <- start(geneCoords) - 20
      promStart <- start(geneCoords) - 2000
    }

    promoter_gr <- GRanges(seqnames = geneChrom, ranges = IRanges(promStart, promEnd))

    if (regulatorType == "proximal") {
      result_list <- lapply(result_list, function(gr) {
        gr[overlapsAny(gr, promoter_gr, type = "within")]
      })
    } else {
      result_list <- lapply(result_list, function(gr) {
        gr[!overlapsAny(gr, promoter_gr, type = "within")]
      })
    }

    result_list <- result_list[vapply(result_list, length, integer(1)) > 0]
  }

  if (regulatorType == "proximal") {
    object@misc$motifs.inProximalPeaks[[key]] <- result_list
  } else {
    object@misc$motifs.inDistalPeaks[[key]] <- result_list
  }

  object
}

# ------------------------------------------------------------
# sequential cache preparation
# ------------------------------------------------------------

prepare_motif_overlap_cache <- function(
    object,
    gene_peak_links = NULL,
    assay = "peaks",
    test.in_promoter = FALSE,
    verbose = TRUE
) {
  if (is.null(gene_peak_links)) {
    gene_peak_links <- object@misc$peak_stats
  }

  helper.df <- unique(gene_peak_links[, c("gene", "cluster"), drop = FALSE])

  if (verbose) {
    message("Preparing motif overlap cache sequentially ...")
    message("n gene/cluster pairs: ", nrow(helper.df))
  }

  for (i in seq_len(nrow(helper.df))) {
    gene <- helper.df$gene[[i]]
    cluster <- helper.df$cluster[[i]]

    if (verbose && (i %% 50 == 0 || i == 1 || i == nrow(helper.df))) {
      message("[cache] ", i, "/", nrow(helper.df), " : ", gene, " | ", cluster)
    }

    object <- compute_motifs_to_gene(
      object = object,
      gene = gene,
      cluster = cluster,
      assay = assay,
      regulatorType = "proximal",
      gene_peak_links = gene_peak_links,
      overwrite = FALSE,
      test.in_promoter = test.in_promoter
    )

    object <- compute_motifs_to_gene(
      object = object,
      gene = gene,
      cluster = cluster,
      assay = assay,
      regulatorType = "distal",
      gene_peak_links = gene_peak_links,
      overwrite = FALSE,
      test.in_promoter = test.in_promoter
    )
  }

  object
}

# ------------------------------------------------------------
# pure worker: returns one data.frame for one gene/cluster
# ------------------------------------------------------------

motif_enrichment_per_gene_df <- function(
    gene,
    cluster,
    gene_peak_links,
    motif_data,
    peak_meta,
    motifs_in_proximal,
    motifs_in_distal,
    motif_granges_cache,
    bg_size = 30
) {
  key <- .make_gene_cluster_key(gene, cluster)

  result_list.proximal <- motifs_in_proximal[[key]]
  result_list.distal   <- motifs_in_distal[[key]]

  if (is.null(result_list.proximal)) result_list.proximal <- list()
  if (is.null(result_list.distal))   result_list.distal   <- list()

  col_names <- colnames(motif_data)

  proximal.peaks_linked_to_gene <- gene_peak_links$peak[
    gene_peak_links$gene == gene &
    gene_peak_links$cluster == cluster &
    gene_peak_links$regulatorType == "proximal"
  ]

  distal.peaks_linked_to_gene <- gene_peak_links$peak[
    gene_peak_links$gene == gene &
    gene_peak_links$cluster == cluster &
    gene_peak_links$regulatorType == "distal"
  ]

  proximal.motif_count <- .count_motifs_in_peaks(motif_data, proximal.peaks_linked_to_gene, col_names)
  distal.motif_count   <- .count_motifs_in_peaks(motif_data, distal.peaks_linked_to_gene, col_names)

  new_motif_counts.proximal <- vapply(result_list.proximal, length, integer(1))
  new_motif_counts.distal   <- vapply(result_list.distal,   length, integer(1))

  if (length(new_motif_counts.proximal) > 0) {
    proximal.motif_count[names(new_motif_counts.proximal)] <- new_motif_counts.proximal
  }
  if (length(new_motif_counts.distal) > 0) {
    distal.motif_count[names(new_motif_counts.distal)] <- new_motif_counts.distal
  }

  calculate_bg_mat <- function(peaks_linked_to_gene) {
    if (length(peaks_linked_to_gene) == 0) {
      return(NULL)
    }

    bg_mat <- matrix(
      0,
      nrow = bg_size,
      ncol = length(col_names),
      dimnames = list(NULL, col_names)
    )

    for (i in seq_len(bg_size)) {
      sampled_peaks_for_gene <- vapply(peaks_linked_to_gene, function(peak) {
        if (!peak %in% rownames(peak_meta)) {
          return(NA_character_)
        }

        cl <- peak_meta[peak, "cluster"]
        comparable_peaks <- rownames(peak_meta)[peak_meta$cluster == cl]

        if (length(comparable_peaks) == 0) {
          return(NA_character_)
        }

        sample(comparable_peaks, 1)
      }, character(1))

      sampled_peaks_for_gene <- sampled_peaks_for_gene[!is.na(sampled_peaks_for_gene)]

      if (length(sampled_peaks_for_gene) == 0) {
        next
      }

      bg_motif_count <- .count_motifs_in_peaks(motif_data, sampled_peaks_for_gene, col_names)
      motifs_for_bg_gene <- names(bg_motif_count)[bg_motif_count > 0]

      if (length(motifs_for_bg_gene) > 0) {
        peaks_GRanges <- .peaks_to_granges(sampled_peaks_for_gene)

        result_list <- lapply(motif_granges_cache[motifs_for_bg_gene], function(gr) {
          subsetByOverlaps(gr, peaks_GRanges)
        })

        new_motif_counts <- vapply(result_list, length, integer(1))
        bg_motif_count[names(new_motif_counts)] <- new_motif_counts
      }

      bg_mat[i, ] <- bg_motif_count
    }

    bg_mat
  }

  compute_stats <- function(fg_count, peaks_linked_to_gene) {
    bg_count <- setNames(rep(0, length(col_names)), col_names)
    log2FC   <- setNames(rep(NA_real_, length(col_names)), col_names)
    t_stat   <- setNames(rep(NA_real_, length(col_names)), col_names)
    p_value  <- setNames(rep(1, length(col_names)), col_names)
    p_adjust <- setNames(rep(1, length(col_names)), col_names)

    if (length(peaks_linked_to_gene) > 0) {
      bg_mat <- calculate_bg_mat(peaks_linked_to_gene)

      if (!is.null(bg_mat) && nrow(bg_mat) > 0) {
        bg_count <- colMeans(bg_mat)

        tt <- vapply(seq_len(ncol(bg_mat)), function(i) {
          .safe_t_test_less(bg_mat[, i], fg_count[[i]])
        }, numeric(2))

        t_stat  <- tt["statistic", ]
        p_value <- tt["p.value", ]
        p_adjust <- pmin(1, p_value * length(col_names))

        suppressWarnings({
          log2FC <- log2(fg_count / bg_count)
        })
      }
    }

    list(
      background_count = bg_count,
      log2FC = log2FC,
      t_stat = t_stat,
      p_value = p_value,
      p_adjust = p_adjust
    )
  }

  prox_stats <- compute_stats(proximal.motif_count, proximal.peaks_linked_to_gene)
  dist_stats <- compute_stats(distal.motif_count, distal.peaks_linked_to_gene)

  data.frame(
    gene = gene,
    cluster = cluster,
    motif = substr(col_names, 1, 8),

    proximal.motif_count = as.numeric(proximal.motif_count),
    proximal.background_count = as.numeric(prox_stats$background_count),
    log2FC.proximal = as.numeric(prox_stats$log2FC),
    t_stat.proximal = as.numeric(prox_stats$t_stat),
    p_value.proximal = as.numeric(prox_stats$p_value),
    p_adjust.proximal = as.numeric(prox_stats$p_adjust),

    distal.motif_count = as.numeric(distal.motif_count),
    distal.background_count = as.numeric(dist_stats$background_count),
    log2FC.distal = as.numeric(dist_stats$log2FC),
    t_stat.distal = as.numeric(dist_stats$t_stat),
    p_value.distal = as.numeric(dist_stats$p_value),
    p_adjust.distal = as.numeric(dist_stats$p_adjust),

    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------
# sequential compatibility wrapper
# keeps old style: returns object with misc$motif_enrichment
# ------------------------------------------------------------

motif_enrichment_per_gene <- function(
    object, gene, cluster, gene_peak_links = NULL,
    assay = "peaks", bg_size = 30,
    test.in_promoter = FALSE
) {
  if (is.null(gene_peak_links)) {
    gene_peak_links <- object@misc$peak_stats
  }

  object <- prepare_motif_overlap_cache(
    object = object,
    gene_peak_links = unique(gene_peak_links[
      gene_peak_links$gene == gene & gene_peak_links$cluster == cluster,
      c("gene", "cluster"),
      drop = FALSE
    ]),
    assay = assay,
    test.in_promoter = test.in_promoter,
    verbose = FALSE
  )

  motif_data <- object[[assay]]@motifs@data
  peak_meta  <- object[[assay]]@meta.features

  motif_granges_cache <- get_GRanges_for_motifs(object, colnames(motif_data))

  res_df <- motif_enrichment_per_gene_df(
    gene = gene,
    cluster = cluster,
    gene_peak_links = gene_peak_links,
    motif_data = motif_data,
    peak_meta = peak_meta,
    motifs_in_proximal = object@misc$motifs.inProximalPeaks,
    motifs_in_distal   = object@misc$motifs.inDistalPeaks,
    motif_granges_cache = motif_granges_cache,
    bg_size = bg_size
  )

  if (is.null(object@misc$motif_enrichment)) {
    object@misc$motif_enrichment <- res_df
  } else {
    object@misc$motif_enrichment <- .rbind_fill_df(object@misc$motif_enrichment, res_df)
  }

  object
}

# ------------------------------------------------------------
# promoter length helper
# ------------------------------------------------------------

.promoter_length_for_gene <- function(gene, peak_stats_filtered, links_df) {
  prox.peak <- peak_stats_filtered$peak[
    peak_stats_filtered$gene == gene &
    peak_stats_filtered$regulatorType == "proximal"
  ]

  prox.peak.dist <- links_df[
    links_df$gene == gene & links_df$peak %in% prox.peak,
  ]$distancePeakStartTSS

  prox.peak.dist <- as.numeric(prox.peak.dist)

  if (length(prox.peak.dist) > 0 && min(prox.peak.dist) < -2000) {
    message("extended promoter")
    return(-min(prox.peak.dist))
  } else {
    return(2000)
  }
}

# ------------------------------------------------------------
# PARALLEL calculate_enrichments
# ------------------------------------------------------------

calculate_enrichments <- function(
    object,
    gene_peak_links = NULL,
    add_to_obj = TRUE,
    bg_size = 50,
    assay = "peaks",
    workers = NULL,
    future_strategy = NULL,
    test.in_promoter = FALSE,
    verbose = TRUE
) {
  # Args:
  #   object: Seurat object
  #   gene_peak_links: data.frame; default object@misc$peak_stats
  #   add_to_obj: if TRUE, writes object@misc$motif_enrichment and returns object
  #   bg_size: number of random bg samples
  #   assay: ChromatinAssay
  #   workers: number of workers; if NULL choose sensible default
  #   future_strategy: optional, e.g. "multisession" or "multicore"
  #   test.in_promoter: forwarded to cache building
  #
  # Returns:
  #   object or motif_enrichment data.frame

  if (is.null(gene_peak_links)) {
    message("using peak stats from object")
    gene_peak_links <- object@misc$peak_stats
  }

  helper.df <- unique(gene_peak_links[, c("gene", "cluster"), drop = FALSE])

  if (verbose) {
    message("counting motifs in linked peaks")
    message("n gene/cluster pairs: ", nrow(helper.df))
  }

  # 1) build motif overlap cache sequentially
  object <- prepare_motif_overlap_cache(
    object = object,
    gene_peak_links = gene_peak_links,
    assay = assay,
    test.in_promoter = test.in_promoter,
    verbose = verbose
  )

  # 2) extract only what workers need
  motif_data <- object[[assay]]@motifs@data
  peak_meta  <- object[[assay]]@meta.features
  motifs_in_proximal <- object@misc$motifs.inProximalPeaks
  motifs_in_distal   <- object@misc$motifs.inDistalPeaks

  if (verbose) {
    message("precomputing motif GRanges cache")
  }
  motif_granges_cache <- get_GRanges_for_motifs(object, colnames(motif_data))

  # 3) set future plan
  if (is.null(workers)) {
    workers <- max(1, min(4, future::availableCores()))
  }

  if (is.null(future_strategy)) {
    future_strategy <- if (.Platform$OS.type == "windows") "multisession" else "multicore"
  }

  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)

  if (future_strategy == "multicore") {
    future::plan(future::multicore, workers = workers)
  } else if (future_strategy == "multisession") {
    future::plan(future::multisession, workers = workers)
  } else {
    stop("future_strategy must be 'multicore' or 'multisession'")
  }

  if (verbose) {
    message("running per-gene motif enrichment in parallel with ", workers, " worker(s)")
  }

    # optional: raise the cap temporarily after shrinking FUN
  oopts <- options(future.globals.maxSize = 2 * 1024^3)  # 2 GiB
  on.exit(options(oopts), add = TRUE)

  res_list <- future.apply::future_lapply(
    X = seq_len(nrow(helper.df)),
    FUN = .motif_worker,
    helper.df = helper.df,
    gene_peak_links = gene_peak_links,
    motif_data = motif_data,
    peak_meta = peak_meta,
    motifs_in_proximal = motifs_in_proximal,
    motifs_in_distal = motifs_in_distal,
    motif_granges_cache = motif_granges_cache,
    bg_size = bg_size,
    future.seed = TRUE,
    future.globals = list(
      .motif_worker = .motif_worker,
      motif_enrichment_per_gene_df = motif_enrichment_per_gene_df,
      .make_gene_cluster_key = .make_gene_cluster_key,
      .peaks_to_granges = .peaks_to_granges,
      .count_motifs_in_peaks = .count_motifs_in_peaks,
      .safe_t_test_less = .safe_t_test_less
    ),
    future.packages = c("GenomicRanges", "IRanges")
  )

  motif_enrichment <- do.call(rbind, res_list)

  # 4) promoter motif counts
  if (verbose) {
    message("counting motifs in promoter region")
  }

  peak_stats_filtered <- object@misc$peak_stats.filtered
  links_df <- Links(object)

  genes <- unique(motif_enrichment$gene)
  prom.length_of.gene <- lapply(genes, function(gene) {
    .promoter_length_for_gene(gene, peak_stats_filtered = peak_stats_filtered, links_df = links_df)
  })
  names(prom.length_of.gene) <- genes

  promoter_counts <- vapply(seq_len(nrow(motif_enrichment)), function(i) {
    gene <- motif_enrichment$gene[[i]]
    motif <- motif_enrichment$motif[[i]]
    prom.length <- prom.length_of.gene[[gene]]
    get_tf_bindingsites_in_region(gene.name = gene, motif.id = motif, prom.length = prom.length)
  }, numeric(1))

  motif_enrichment$promoter.motif_count <- promoter_counts

  if (add_to_obj) {
    object@misc$motif_enrichment <- motif_enrichment
    return(object)
  } else {
    return(motif_enrichment)
  }
}

getMotifEnrichment <- function(object, gene, cluster, motif = NULL) {
  if (is.null(object@misc$motif_enrichment)) {
    object <- motif_enrichment_per_gene(object = object, gene = gene, cluster = cluster)
  }

  if (!is.null(motif)) {
    return(
      object@misc$motif_enrichment[
        object@misc$motif_enrichment[, "gene"] == gene &
        object@misc$motif_enrichment[, "motif"] == motif &
        object@misc$motif_enrichment[, "cluster"] == cluster,
        ,
        drop = FALSE
      ]
    )
  }

  object@misc$motif_enrichment[
    object@misc$motif_enrichment[, "gene"] == gene &
    object@misc$motif_enrichment[, "cluster"] == cluster,
    ,
    drop = FALSE
  ]
}

getMotifStats <- function(object, gene, cluster, motif = NULL, gene_peak_links = NULL) {
  if (!is.null(object@misc$motif_stats)) {
    if (!is.null(motif)) {
      motif_stats_for_seed <- object@misc$motif_stats[
        object@misc$motif_stats[, "gene"] == gene &
        object@misc$motif_stats[, "motif"] == motif &
        object@misc$motif_stats[, "cluster"] == cluster,
        ,
        drop = FALSE
      ]
      return(motif_stats_for_seed)
    }

    motif_stats_for_seed <- object@misc$motif_stats[
      object@misc$motif_stats[, "gene"] == gene &
      object@misc$motif_stats[, "cluster"] == cluster,
      ,
      drop = FALSE
    ]

    if (nrow(motif_stats_for_seed) > 0) {
      return(motif_stats_for_seed)
    }

    return(NULL)
  }

  object <- motif_enrichment_per_gene(
    object = object,
    gene = gene,
    cluster = cluster,
    gene_peak_links = gene_peak_links
  )

  object@misc$motif_enrichment[
    object@misc$motif_enrichment[, "gene"] == gene &
    object@misc$motif_enrichment[, "cluster"] == cluster,
    ,
    drop = FALSE
  ]
}