@export
plot_module_scores <- function(
  object,
  seed_genes,
  module_name = "AVN_module",
  assay_name = "RNA",
  cluster_name = "cell_type",
  test_method = c("wilcox", "t.test"),
  alternative = c("greater", "two.sided", "less"),
  p_adjust_method = "BH",
  label_type = c("stars", "p.adj"),
  show_ns = TRUE
) {
  test_method <- match.arg(test_method)
  alternative <- match.arg(alternative)
  label_type <- match.arg(label_type)

  DefaultAssay(object) <- assay_name
  score_col <- paste0(module_name, "1")

  # Keep only genes present in this object
  seed_genes_present <- intersect(seed_genes, rownames(object))
  missing_genes <- setdiff(seed_genes, seed_genes_present)

  if (length(seed_genes_present) == 0) {
    stop("None of the seed genes are present in the object.")
  }

  cat("Genes found:", length(seed_genes_present), "\n")
  if (length(missing_genes) > 0) {
    cat("Missing genes:", paste(missing_genes, collapse = ", "), "\n")
  }

  Idents(object) <- cluster_name

  # Add module score
  object <- AddModuleScore(
    object = object,
    features = list(seed_genes_present),
    name = module_name,
    assay = assay_name,
    seed = 123
  )

  # Pull data for plotting
  plot_df <- FetchData(object, vars = c(cluster_name, score_col))
  colnames(plot_df) <- c("cluster", "score")

  # Order clusters by median score
  celltype_order <- plot_df %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarise(med = median(score, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(med)) %>%
    dplyr::pull(cluster)

  plot_df$cluster <- factor(plot_df$cluster, levels = celltype_order)

  # Per-cluster enrichment test: cluster vs all other cells
  stats_df <- lapply(levels(plot_df$cluster), function(cl) {
    x <- plot_df$score[plot_df$cluster == cl]
    y <- plot_df$score[plot_df$cluster != cl]

    test_res <- switch(
      test_method,
      wilcox = stats::wilcox.test(x, y, alternative = alternative, exact = FALSE),
      `t.test` = stats::t.test(x, y, alternative = alternative)
    )

    data.frame(
      cluster = cl,
      n_cluster = length(x),
      n_other = length(y),
      mean_cluster = mean(x, na.rm = TRUE),
      mean_other = mean(y, na.rm = TRUE),
      median_cluster = median(x, na.rm = TRUE),
      median_other = median(y, na.rm = TRUE),
      effect_mean = mean(x, na.rm = TRUE) - mean(y, na.rm = TRUE),
      effect_median = median(x, na.rm = TRUE) - median(y, na.rm = TRUE),
      p_val = test_res$p.value,
      stringsAsFactors = FALSE
    )
  }) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(
      p_adj = stats::p.adjust(p_val, method = p_adjust_method)
    )

  # Significance labels
  p_to_stars <- function(p) {
    dplyr::case_when(
      is.na(p)        ~ NA_character_,
      p <= 1e-4       ~ "****",
      p <= 1e-3       ~ "***",
      p <= 1e-2       ~ "**",
      p <= 5e-2       ~ "*",
      TRUE            ~ "ns"
    )
  }

  stats_df <- stats_df %>%
    dplyr::mutate(
      label = if (label_type == "stars") {
        p_to_stars(p_adj)
      } else {
        paste0("adj.P=", signif(p_adj, 2))
      }
    )

  if (!show_ns && label_type == "stars") {
    stats_df$label[stats_df$label == "ns"] <- ""
  }

  # Y-positions for annotation
  y_range <- diff(range(plot_df$score, na.rm = TRUE))
  if (y_range == 0) y_range <- 1

  annot_df <- plot_df %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarise(y_pos = max(score, na.rm = TRUE) + 0.06 * y_range, .groups = "drop") %>%
    dplyr::left_join(stats_df, by = "cluster")

  # Plot
  p <- ggplot(plot_df, aes(x = cluster, y = score, fill = cluster)) +
    geom_violin(trim = FALSE, scale = "width", color = "grey30") +
    geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.5, color = "white") +
    geom_text(
      data = annot_df,
      aes(x = cluster, y = y_pos, label = label),
      inherit.aes = FALSE,
      vjust = 0,
      size = 4
    ) +
    theme_classic() +
    labs(
      x = cluster_name,
      y = score_col,
      subtitle = paste0(
        "Per-cluster enrichment test: ",
        test_method,
        " (", alternative, "), p-adjust = ", p_adjust_method
      )
    ) +
    scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.15))) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )

  return(list(
    plot = p,
    object = object,
    stats = stats_df,
    plot_data = plot_df
  ))
}