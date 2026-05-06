#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(jsonlite)
})

in_dir   <- "/workspaces/mopitas-mapra/backend/dgea/export_for_seurat"
out_json <- "/workspaces/mopitas-mapra/backend/dgea/results/dgea_results.json"

wanted_obs_cols <- c("cell_type", "leiden", "niche_cluster")

# Compare all group pairs within each obs column
do_all_vs_all <- TRUE
do_one_vs_all <- TRUE

# FindMarkers Parameter
min_cells <- 3
top_n <- 200
min_pct <- 0.10
logfc_threshold <- 0.25


# Heatmap-Parameter
heatmap_n_top <- 10
heatmap_padj_cutoff <- 0.05
heatmap_logfc_cutoff <- 0.5
heatmap_slot <- "data"   # normalized data for heatmap
# optional: filter cells to "in_tissue" == 1
filter_in_tissue <- TRUE

wanted_marker_cols <- c("avg_log2FC", "avg_logFC", "p_val", "p_val_adj", "pct1", "pct2")

# =========================
# Helper
# =========================
safe_id <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

# Convert df (gene + cols) -> TableData {col: {gene: value}}
df_to_tabledata <- function(df, index_col = "gene") {
  out <- list()
  cols <- setdiff(colnames(df), index_col)
  for (c in cols) out[[c]] <- list()

  for (i in seq_len(nrow(df))) {
    idx <- as.character(df[[index_col]][i])
    for (c in cols) {
      v <- df[[c]][i]
      if (!is.na(v) && !is.null(v)) {
        out[[c]][[idx]] <- v
      }
    }
  }
  out
}

# robust numeric conversion
to_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

build_heatmap_context <- function(
  seu,
  obs_col,
  group1,
  group2,
  table_df,
  assay = "RNA",
  slot = "data",
  n_top = 10,
  padj_cutoff = 0.05,
  logfc_cutoff = 0.5
) {
  fc_col <- NULL
  if ("avg_log2FC" %in% colnames(table_df)) {
    fc_col <- "avg_log2FC"
  } else if ("avg_logFC" %in% colnames(table_df)) {
    fc_col <- "avg_logFC"
  } else {
    return(list(groups = list(), rows = list()))
  }

  df <- table_df

  # Filter by adjusted p-value if available
  if ("p_val_adj" %in% colnames(df)) {
    df <- df[is.na(df$p_val_adj) | df$p_val_adj < padj_cutoff, , drop = FALSE]
  }

  # Filter by log fold change
  df <- df[abs(df[[fc_col]]) >= logfc_cutoff, , drop = FALSE]

  if (nrow(df) == 0) {
    return(list(groups = list(), rows = list()))
  }

  # # Top genes upregulated in group1 and group2
  pos <- df[df[[fc_col]] > 0, , drop = FALSE]
  neg <- df[df[[fc_col]] < 0, , drop = FALSE]

  pos <- pos[order(-pos[[fc_col]]), , drop = FALSE]
  neg <- neg[order(neg[[fc_col]]), , drop = FALSE]

  top_pos <- head(pos$gene, n_top)
  top_neg <- head(neg$gene, n_top)

  genes_use <- unique(c(top_pos, top_neg))

  if (length(genes_use) == 0) {
    return(list(groups = list(), rows = list()))
  }

  expr_mat <- GetAssayData(seu, assay = assay, slot = slot)
  group_values <- as.character(seu[[obs_col, drop = TRUE]])

  all_groups <- sort(unique(group_values))
  all_groups <- all_groups[!is.na(all_groups)]

  # first two groups are the ones compared, then the rest (for consistent column order in heatmap)
  other_groups <- setdiff(all_groups, c(group1, group2))
  groups_ordered <- c(group1, group2, other_groups)

  mat_raw <- sapply(groups_ordered, function(g) {
    cells <- colnames(seu)[group_values == g]
    if (length(cells) == 0) {
      rep(NA_real_, length(genes_use))
    } else {
      Matrix::rowMeans(expr_mat[genes_use, cells, drop = FALSE])
    }
  })

  mat_raw <- as.matrix(mat_raw)
  rownames(mat_raw) <- genes_use
  colnames(mat_raw) <- groups_ordered

  # row scaling
  mat_scaled <- t(scale(t(mat_raw)))
  mat_scaled[is.na(mat_scaled)] <- 0

  rows <- lapply(seq_len(nrow(mat_raw)), function(i) {
    list(
      gene = rownames(mat_raw)[i],
      scaled = unname(as.numeric(mat_scaled[i, ])),
      raw = unname(as.numeric(mat_raw[i, ]))
    )
  })

  list(
    groups = as.list(colnames(mat_raw)),
    rows = rows
  )
}

build_heatmap_context_vs_all <- function(
  seu,
  obs_col,
  group1,
  table_df,
  assay = "RNA",
  slot = "data",
  n_top = 10,
  padj_cutoff = 0.05,
  logfc_cutoff = 0.5
) {
  fc_col <- NULL
  if ("avg_log2FC" %in% colnames(table_df)) {
    fc_col <- "avg_log2FC"
  } else if ("avg_logFC" %in% colnames(table_df)) {
    fc_col <- "avg_logFC"
  } else {
    return(list(groups = list(), rows = list()))
  }

  df <- table_df

  if ("p_val_adj" %in% colnames(df)) {
    df <- df[is.na(df$p_val_adj) | df$p_val_adj < padj_cutoff, , drop = FALSE]
  }

  df <- df[abs(df[[fc_col]]) >= logfc_cutoff, , drop = FALSE]

  if (nrow(df) == 0) {
    return(list(groups = list(), rows = list()))
  }

  pos <- df[df[[fc_col]] > 0, , drop = FALSE]
  neg <- df[df[[fc_col]] < 0, , drop = FALSE]

  pos <- pos[order(-pos[[fc_col]]), , drop = FALSE]
  neg <- neg[order(neg[[fc_col]]), , drop = FALSE]

  top_pos <- head(pos$gene, n_top)
  top_neg <- head(neg$gene, n_top)

  genes_use <- unique(c(top_pos, top_neg))

  if (length(genes_use) == 0) {
    return(list(groups = list(), rows = list()))
  }

  expr_mat <- GetAssayData(seu, assay = assay, slot = slot)
  group_values <- as.character(seu[[obs_col, drop = TRUE]])

  all_groups <- sort(unique(group_values))
  all_groups <- all_groups[!is.na(all_groups)]

  other_groups <- setdiff(all_groups, group1)
  groups_ordered <- c(group1, other_groups)

  mat_raw <- sapply(groups_ordered, function(g) {
    cells <- colnames(seu)[group_values == g]
    if (length(cells) == 0) {
      rep(NA_real_, length(genes_use))
    } else {
      Matrix::rowMeans(expr_mat[genes_use, cells, drop = FALSE])
    }
  })

  mat_raw <- as.matrix(mat_raw)
  rownames(mat_raw) <- genes_use
  colnames(mat_raw) <- groups_ordered

  mat_scaled <- t(scale(t(mat_raw)))
  mat_scaled[is.na(mat_scaled)] <- 0

  rows <- lapply(seq_len(nrow(mat_raw)), function(i) {
    list(
      gene = rownames(mat_raw)[i],
      scaled = unname(as.numeric(mat_scaled[i, ])),
      raw = unname(as.numeric(mat_raw[i, ]))
    )
  })

  list(
    groups = as.list(colnames(mat_raw)),
    rows = rows
  )
}

compute_one_vs_all_tabledata <- function(seu, obs_col, g1) {
  Idents(seu) <- obs_col
  sizes <- table(Idents(seu))
  n1 <- as.integer(sizes[[g1]]); if (is.na(n1)) n1 <- 0L
  n2 <- as.integer(sum(sizes)) - n1

  cmp_id <- paste0(safe_id(g1), "__vs__all")
  g2_label <- "all"

  if (n1 < min_cells || n2 < min_cells) {
    return(list(
      id = cmp_id,
      name = paste(g1, "vs", g2_label),
      group1 = g1,
      group2 = g2_label,
      n1 = n1,
      n2 = n2,
      skipped = TRUE,
      skip_reason = sprintf("min_cells=%d not met", min_cells),
      table = list(),
      heatmap_context = list(groups = list(), rows = list())
    ))
  }

  others <- setdiff(levels(Idents(seu)), g1)

  res <- tryCatch({
    FindMarkers(
      seu,
      ident.1 = g1,
      ident.2 = others,
      test.use = "wilcox",
      min.pct = min_pct,
      logfc.threshold = logfc_threshold
    )
  }, error = function(e) {
    structure(NULL, error_msg = e$message)
  })

  if (is.null(res)) {
    return(list(
      id = cmp_id,
      name = paste(g1, "vs", g2_label),
      group1 = g1,
      group2 = g2_label,
      n1 = n1,
      n2 = n2,
      skipped = TRUE,
      skip_reason = paste0("FindMarkers error: ", attr(res, "error_msg")),
      table = list(),
      heatmap_context = list(groups = list(), rows = list())
    ))
  }

  res$gene <- rownames(res)

  if ("p_val_adj" %in% colnames(res)) {
    res <- res[order(res$p_val_adj, res$p_val), , drop = FALSE]
  } else if ("p_val" %in% colnames(res)) {
    res <- res[order(res$p_val), , drop = FALSE]
  }
  if (nrow(res) > top_n) res <- res[1:top_n, , drop = FALSE]

  if ("pct.1" %in% colnames(res)) res$pct1 <- res[["pct.1"]]
  if ("pct.2" %in% colnames(res)) res$pct2 <- res[["pct.2"]]

  keep <- intersect(wanted_marker_cols, colnames(res))
  table_df <- res[, c("gene", keep), drop = FALSE]

  for (c in keep) table_df[[c]] <- to_num(table_df[[c]])

  heatmap_context <- build_heatmap_context_vs_all(
    seu = seu,
    obs_col = obs_col,
    group1 = g1,
    table_df = table_df,
    assay = DefaultAssay(seu),
    slot = heatmap_slot,
    n_top = heatmap_n_top,
    padj_cutoff = heatmap_padj_cutoff,
    logfc_cutoff = heatmap_logfc_cutoff
  )

  list(
    id = cmp_id,
    name = paste(g1, "vs", g2_label),
    group1 = g1,
    group2 = g2_label,
    n1 = n1,
    n2 = n2,
    skipped = FALSE,
    table = df_to_tabledata(table_df, index_col = "gene"),
    heatmap_context = heatmap_context
  )
}

compute_one_comparison_tabledata <- function(seu, obs_col, g1, g2) {
  Idents(seu) <- obs_col
  sizes <- table(Idents(seu))
  n1 <- as.integer(sizes[[g1]]); if (is.na(n1)) n1 <- 0L
  n2 <- as.integer(sizes[[g2]]); if (is.na(n2)) n2 <- 0L

  cmp_id <- paste0(safe_id(g1), "__vs__", safe_id(g2))

  # skip if not enough cells in either group (FindMarkers requirement)
  if (n1 < min_cells || n2 < min_cells) {
    return(list(
      id = cmp_id,
      name = paste(g1, "vs", g2),
      group1 = g1,
      group2 = g2,
      n1 = n1,
      n2 = n2,
      skipped = TRUE,
      skip_reason = sprintf("min_cells=%d not met", min_cells),
      table = list(),
      heatmap_context = list(groups = list(), rows = list())
    ))
  }

  res <- tryCatch({
    FindMarkers(
      seu,
      ident.1 = g1,
      ident.2 = g2,
      test.use = "wilcox",
      min.pct = min_pct,
      logfc.threshold = logfc_threshold
    )
  }, error = function(e) {
    structure(NULL, error_msg = e$message)
  })

  if (is.null(res)) {
    return(list(
      id = cmp_id,
      name = paste(g1, "vs", g2),
      group1 = g1,
      group2 = g2,
      n1 = n1,
      n2 = n2,
      skipped = TRUE,
      skip_reason = paste0("FindMarkers error: ", attr(res, "error_msg")),
      table = list(),
      heatmap_context = list(groups = list(), rows = list())
    ))
  }

  res$gene <- rownames(res)

  # sort + top_n
  if ("p_val_adj" %in% colnames(res)) {
    res <- res[order(res$p_val_adj, res$p_val), , drop = FALSE]
  } else if ("p_val" %in% colnames(res)) {
    res <- res[order(res$p_val), , drop = FALSE]
  }
  if (nrow(res) > top_n) res <- res[1:top_n, , drop = FALSE]

  # pct.1/pct.2 -> pct1/pct2
  if ("pct.1" %in% colnames(res)) res$pct1 <- res[["pct.1"]]
  if ("pct.2" %in% colnames(res)) res$pct2 <- res[["pct.2"]]

  keep <- intersect(wanted_marker_cols, colnames(res))
  table_df <- res[, c("gene", keep), drop = FALSE]

  # Ensure numeric columns are stored as numeric values
  for (c in keep) table_df[[c]] <- to_num(table_df[[c]])

  heatmap_context <- build_heatmap_context(
    seu = seu,
    obs_col = obs_col,
    group1 = g1,
    group2 = g2,
    table_df = table_df,
    assay = DefaultAssay(seu),
    slot = heatmap_slot,
    n_top = heatmap_n_top,
    padj_cutoff = heatmap_padj_cutoff,
    logfc_cutoff = heatmap_logfc_cutoff
  )

  list(
    id = cmp_id,
    name = paste(g1, "vs", g2),
    group1 = g1,
    group2 = g2,
    n1 = n1,
    n2 = n2,
    skipped = FALSE,
    table = df_to_tabledata(table_df, index_col = "gene"),
    heatmap_context = heatmap_context
  )
}

# Aggregate results across all comparisons to get the best p-value for each gene, for an overall "top markers" view
build_aggregates_tabledata <- function(comparisons) {
  best_rows <- new.env(parent = emptyenv())

  for (cmp in comparisons) {
    if (isTRUE(cmp$skipped)) next
    td <- cmp$table
    if (length(td) == 0) next

    # genes from first column map
    first_col <- names(td)[1]
    genes <- names(td[[first_col]])

    for (g in genes) {
      row <- list(gene = g)
      for (col in names(td)) row[[col]] <- td[[col]][[g]]

      padj <- if (!is.null(row$p_val_adj)) to_num(row$p_val_adj) else NA_real_

      if (!exists(g, envir = best_rows, inherits = FALSE)) {
        assign(g, row, envir = best_rows)
      } else {
        current <- get(g, envir = best_rows, inherits = FALSE)
        current_padj <- if (!is.null(current$p_val_adj)) to_num(current$p_val_adj) else NA_real_
        if (!is.na(padj) && (is.na(current_padj) || padj < current_padj)) {
          assign(g, row, envir = best_rows)
        }
      }
    }
  }

  genes <- ls(best_rows)
  if (length(genes) == 0) return(list())

  df <- do.call(rbind, lapply(genes, function(g) {
    as.data.frame(get(g, envir = best_rows, inherits = FALSE), stringsAsFactors = FALSE)
  }))

  # ensure numeric cols
  for (c in setdiff(colnames(df), "gene")) df[[c]] <- to_num(df[[c]])

  df_to_tabledata(df, index_col = "gene")
}

run_for_obs_col <- function(seu, obs_col) {
  Idents(seu) <- obs_col
  lvls <- levels(Idents(seu))
  sizes <- sort(table(Idents(seu)), decreasing = TRUE)

  pairs <- if (do_all_vs_all) combn(lvls, 2, simplify = FALSE) else list()
  pair_comps <- lapply(pairs, function(p) {
    compute_one_comparison_tabledata(seu, obs_col, p[[1]], p[[2]])
  })

  one_vs_all_comps <- if (do_one_vs_all) {
    lapply(lvls, function(g) compute_one_vs_all_tabledata(seu, obs_col, g))
  } else {
    list()
  }

  comps_list <- c(pair_comps, one_vs_all_comps)

  comparisons <- list()
  for (c in comps_list) comparisons[[c$id]] <- c

  list(
    levels = as.list(lvls),
    group_sizes = setNames(as.list(as.integer(sizes)), names(sizes)),
    aggregates = build_aggregates_tabledata(comps_list),
    comparisons = comparisons
  )
}

mat <- ReadMtx(
  mtx = file.path(in_dir, "matrix.mtx"),
  features = file.path(in_dir, "genes.tsv"),
  cells = file.path(in_dir, "barcodes.tsv"),
  feature.column = 1
)
meta <- read.csv(file.path(in_dir, "metadata.csv"), row.names = 1, check.names = FALSE)

seu <- CreateSeuratObject(counts = mat, meta.data = meta, assay = "RNA")
DefaultAssay(seu) <- "RNA"

if (filter_in_tissue && ("in_tissue" %in% colnames(seu@meta.data))) {
  seu <- subset(seu, subset = in_tissue == 1)
}

seu <- NormalizeData(seu)

# =========================
# Run
# =========================
available_cols <- wanted_obs_cols[wanted_obs_cols %in% colnames(seu@meta.data)]
message("Available obs cols: ", paste(available_cols, collapse = ", "))

dgea <- list()
for (col in available_cols) {
  message("Running DGEA for obs_col: ", col)
  dgea[[col]] <- run_for_obs_col(seu, col)
}

out <- list(
  meta = list(
    dgea = dgea,
    dgea_meta = list(
      method = "Seurat::FindMarkers",
      obs_cols_requested = as.list(wanted_obs_cols),
      obs_cols_available = as.list(available_cols),
      min_cells = min_cells,
      top_n = top_n,
      min_pct = min_pct,
      logfc_threshold = logfc_threshold,
      filter_in_tissue = filter_in_tissue
    )
  )
)

dir.create(dirname(out_json), recursive = TRUE, showWarnings = FALSE)
write_json(out, out_json, auto_unbox = TRUE, pretty = TRUE)
message("Wrote JSON: ", out_json)
