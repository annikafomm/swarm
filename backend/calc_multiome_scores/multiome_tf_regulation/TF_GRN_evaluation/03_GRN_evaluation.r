library(parallel)
library(GenomicFeatures)
set.seed(1234)
library(repr)
library(motifmatchr)
library(fastmatch)
library(qlcMatrix)

library(ensembldb)
library(EnsDb.Hsapiens.v86)
library(AnnotationFilter)

library(Signac)
library(Seurat)
library(JASPAR2024)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(patchwork)
library(ggplot2)
library(Matrix)
library(zoo)
library(tidyr)

source("../multiome_methods/function_calls.r")
source("../multiome_methods/peak_gene_relations.r")
source("../multiome_methods/binding_site_identification.r")
source("../multiome_methods/binding_site_identification_parallel.r")
source("../multiome_methods/TF_activity.r")
source("../multiome_methods/signac_utils.r")

library(dplyr)
library(ggplot2)
library(hexbin)
library(RColorBrewer)
library(cowplot)
library(gridExtra)
library(ggExtra)
library(pheatmap)
library(purrr)
library(data.table)
library(Rsamtools)
library(stringi)
library(jsonlite)
library(optparse)

# ==============================================================================
# GRN EVALUATION PIPELINE - MAIN FUNCTIONS
# ==============================================================================

create_overview_plots <- function(object, ident, analysis_name, outpath) {
  dir.create(paste0(outpath, "/peak_plots"), recursive = TRUE, showWarnings = FALSE)
  
  # UMAP plot
  p <- DimPlot(object, group.by = ident, label = TRUE) 
  ggsave(paste0(outpath, "/peak_plots/umap.png"), plot = p, width = 10, height = 6)
  
  # Module enrichment plot
  res <- plot_module_scores(
    object = object,
    seed_genes = seed_genes,
    module_name = analysis_name,
    assay_name = "RNA",
    cluster_name = "cell_type",
    test_method = "wilcox",
    alternative = "greater",
    label_type = "stars",
    show_ns = FALSE
  )
  ggsave(paste0(outpath, "/peak_plots/module_enrichment.pdf"), plot = res$plot, width = 8, height = 6)
}

compute_grn_analysis <- function(object, seed_genes, priorGRN, cluster, outpath, analysis_name) {
  # 1. Link peaks to genes
  res_links <- linking_peaks_for_module_cluster(
    object, seed_genes, priorGRN, cluster = cluster
  )
  ggsave(paste0(outpath, "/peak_plots/links_distr.png"), plot = res_links$plot_links, width = 6, height = 6)
  
  # 2. Peak stats
  seeds_for_peaks_stats <- data.frame(
    cluster = cluster,
    gene = res_links$seeds[["gene"]],
    annotation = analysis_name,
    stringsAsFactors = FALSE
  )
  res_peak_stats <- do_peak_stats(object = res_links$object, seeds = seeds_for_peaks_stats, promoter.needed = FALSE)
  
  # 3. Motif stats
  res_local_motif_stats <- do_local_motif_stats(
    object = res_peak_stats,
    seeds = seeds_for_peaks_stats,
    fragpath = NULL,
    peaks_assay = 'peaks',
    background_size = 50,
    parallel = TRUE,
    min_count = 1,
    min_log2FC = 1,
    max_t_stat = 0,
    min_p_adjust = 0.05
  )
  
  # 4. Map motifs and build network
  res_motifs_mapped <- add_motif_TF_mapping(res_local_motif_stats)
  colnames(priorGRN)[colnames(priorGRN) == "gene_name"] <- "gene"
  
  object_w_net <- draw_subNetwork_for_context(
    object = res_motifs_mapped,
    gene_list = seed_genes,
    cluster.specific = TRUE,
    priorGRN = priorGRN,
    cluster = cluster,
    context_name = analysis_name,
    resolution = 20,
    pointsize = 50,
    sd = 5,
    score = -0.5
  )
  
  return(object_w_net)
}

export_results <- function(object, outpath, analysis_name) {
  # Save main object
  saveRDS(object, file.path(outpath, "GRN_object.rds"))
  
  # Export peak stats
  summary_col_labels_short <- c(
    gene = "Gene", cluster = "Cluster", annotation = "Annotation", peak = "Peak",
    regulatorType = "Class", signac.scores = "Link Score", signac.zscores = "Link Z",
    signac.pvalues = "Link P", `t-stat_highly.acc` = "Acc. T-stat",
    `p.value-t.test_highly.acc_BH` = "Acc. FDR", acc_cells_cluster = "Accessible Cells",
    `expr_not_0.given_acc` = "P(expr|acc), cluster", `expr_not_0.given_acc.bg` = "P(expr|acc), bg",
    expr_and_acc_not_0 = "P(expr & acc), cluster", expr_and_acc_not_0.all = "P(expr & acc), all",
    FC.expr_given_acc = "Enrichment, cluster", FC.expr_given_acc.all = "Enrichment, all",
    delta_expr_given_acc.same_peak_bg = "Delta P(expr|acc)",
    promotersLinkedToSeed = "Promoter Peaks", distalPeaksLinkedToSeed = "Distal Peaks",
    pass_type = "Pass Type"
  )
  
  summary_cols <- c(
    "gene", "cluster", "annotation", "peak", "regulatorType", "signac.scores",
    "signac.zscores", "signac.pvalues", "t-stat_highly.acc", "p.value-t.test_highly.acc_BH",
    "acc_cells_cluster", "expr_not_0.given_acc", "expr_not_0.given_acc.bg",
    "expr_and_acc_not_0", "expr_and_acc_not_0.all", "FC.expr_given_acc",
    "FC.expr_given_acc.all", "delta_expr_given_acc.same_peak_bg",
    "promotersLinkedToSeed", "distalPeaksLinkedToSeed", "pass_type"
  )
  
  peak_stats_summary <- object@misc$peak_stats[, summary_cols]
  colnames(peak_stats_summary) <- summary_col_labels_short[summary_cols]
  write.csv(peak_stats_summary, paste0(outpath, "/peak_stats.csv"), row.names = FALSE)
  
  # Export motif stats
  motif_summary_col_labels <- c(
    gene = "Gene", cluster = "Cluster", TF = "TF", motif = "Motif",
    proximal.motif_count = "Prox Motif count", proximal.background_count = "Prox Bg count",
    log2FC.proximal = "Prox Log2FC", p_adjust.proximal = "Prox p-value adj",
    distal.motif_count = "Dist Motif count", distal.background_count = "Dist Bg count",
    log2FC.distal = "Dist Log2FC", p_adjust.distal = "Dist p-value adj",
    promoter.motif_count = "Prom Motif count", footprint_score = "FP Score",
    bg_footprint_mean = "Bg FP Score", footprint.p_value_adj = "FP p-value adj",
    bg_size = "Bg Size", sd.flanks = "Flank sd", bg_sd_mean = "Bg Flank sd",
    left_flank_nonzero_positions = "Left Flank != 0", right_flank_nonzero_positions = "Right Flank != 0"
  )
  
  motif_summary_cols <- c(
    "gene", "cluster", "TF", "motif", "proximal.motif_count", "proximal.background_count",
    "log2FC.proximal", "p_adjust.proximal", "distal.motif_count", "distal.background_count",
    "log2FC.distal", "p_adjust.distal", "promoter.motif_count", "footprint_score",
    "bg_footprint_mean", "footprint.p_value_adj", "bg_size", "sd.flanks", "bg_sd_mean",
    "left_flank_nonzero_positions", "right_flank_nonzero_positions"
  )
  
  motif2TF_df <- data.frame(
    motif = rep(names(object@misc$motif2TF), lengths(object@misc$motif2TF)),
    TF = unlist(object@misc$motif2TF, use.names = FALSE),
    stringsAsFactors = FALSE
  )
  
  tmp <- object@misc$motif_stats %>% left_join(motif2TF_df, by = "motif")
  motif_stats_summary <- tmp[, motif_summary_cols]
  colnames(motif_stats_summary) <- motif_summary_col_labels[motif_summary_cols]
  write.csv(motif_stats_summary, paste0(outpath, "/motif_stats.csv"), row.names = FALSE)
}

export_graph_json <- function(object, outpath, analysis_name) {
  df <- object@misc$context_subNetwork[[analysis_name]]
  
  links <- df %>%
    transmute(
      source = from,
      target = to,
      edge_color = color,
      edge_width = reg_type,
      edge_dash = ifelse(in.prom == 1, "4,4", "0")
    )
  
  from_nodes <- df %>% distinct(name = from, priorTF) %>% mutate(type = "from")
  to_nodes <- df %>% distinct(name = to) %>% mutate(priorTF = NA, type = "to")
  
  nodes <- bind_rows(from_nodes, to_nodes) %>%
    group_by(name) %>%
    summarise(
      priorTF = dplyr::first(na.omit(priorTF)),
      type = ifelse(any(type == "from"), "from", "to"),
      .groups = "drop"
    ) %>%
    mutate(
      node_fill = case_when(
        type == "to" ~ "#cccccc",
        priorTF %in% TRUE ~ "#ff7f0e",
        priorTF %in% FALSE ~ "#1f77b4",
        TRUE ~ "#999999"
      )
    ) %>%
    mutate(id = row_number() - 1)
  
  links <- links %>%
    left_join(nodes %>% select(source = name, source_id = id), by = "source") %>%
    left_join(nodes %>% select(target = name, target_id = id), by = "target") %>%
    transmute(source = source_id, target = target_id, edge_color, edge_width, edge_dash)
  
  graph <- list(nodes = nodes, links = links)
  write_json(graph, paste0(outpath, "/graph.json"), auto_unbox = TRUE, pretty = TRUE)
}

# ==============================================================================
# COMMAND-LINE ARGUMENTS
# ==============================================================================

option_list <- list(
  make_option(c("--analysis_name"), type = "character", default = "collagen_fibral_orga",
              help = "Analysis name for output organization"),
  make_option(c("--gene_set_path"), type = "character", default = "gene_module.csv",
              help = "Path to CSV file with gene module (must have SYMBOL column)"),
  make_option(c("--priorGRN_path"), type = "character", default = "priorGRN.csv",
              help = "Path to prior GRN CSV file"),
  make_option(c("--global_analysis_out_path"), type = "character", default = "global_motif_analysis.rds",
              help = "Path to saved global analysis RDS object"),
  make_option(c("--cluster"), type = "character", default = "Fibroblast",
              help = "Cluster name for analysis"),
  make_option(c("--ident"), type = "character", default = "cell_type",
              help = "Identity variable name"),
  make_option(c("--jobdir_dir"), type = "character", default = "",
              help = "Base job directory for output organization"),
  make_option(c("--initial_upload"), type = "logical", default = TRUE,
              help = "Whether this is initial upload (TRUE) or on-demand (FALSE)")
)

opt <- parse_args(OptionParser(option_list = option_list))

analysis_name <- opt$analysis_name
gene_set_path <- opt$gene_set_path
priorGRN_path <- opt$priorGRN_path
cluster <- opt$cluster
ident <- opt$ident
initial_upload <- opt$initial_upload
global_analysis_out_path <- opt$global_analysis_out_path
jobdir_dir <- opt$jobdir_dir

# Load data
seed_genes <- read.csv(gene_set_path)$SYMBOL
priorGRN <- read.csv(priorGRN_path)
global_analysis_out <- readRDS(global_analysis_out_path)

# Determine output path
if (initial_upload) {
  outpath <- paste0(jobdir_dir, "/multiome/GRN_evaluation/initial_upload/", analysis_name)
} else {
  outpath <- paste0(jobdir_dir, "/multiome/GRN_evaluation/on_demand/", analysis_name)
}

# ==============================================================================
# RUN PIPELINE
# ==============================================================================

create_overview_plots(global_analysis_out, ident, analysis_name, outpath)

object_w_net <- compute_grn_analysis(global_analysis_out, seed_genes, priorGRN, cluster, outpath, analysis_name)

export_results(object_w_net, outpath, analysis_name)

export_graph_json(object_w_net, outpath, analysis_name)




