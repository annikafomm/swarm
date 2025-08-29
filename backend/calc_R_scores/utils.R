calc_aucell_score <- function(expr, modules, n_cores=1) {
  
  # Create AUCell object
  rankings <- AUCell_buildRankings(expr, plotStats=FALSE, verbose=FALSE)
  
  # Calculate AUCell scores for the SPONGE network
  aucell_scores <- AUCell_calcAUC(modules, rankings, aucMaxRank = nrow(expr), nCores = n_cores, verbose=FALSE)
  aucell_scores <- getAUC(aucell_scores)
  aucell_scores <- as.data.frame(aucell_scores)
  
  return(aucell_scores)
}

# TODO: choose parameters wisely
calc_spongeffects_oe <- function(expr, modules, n_cores=1) {
  
  oe.modules <- enrichment_modules(Expr.matrix = expr, modules = modules, bin.size = 10, min.size = 1, max.size = Inf, min.expr = 1, method = "OE", cores=n_cores)
  
  return(oe.modules)
}

calc_spongeffects_gsva <- function(expr, modules, n_cores=1) {
  
  gsvaPar <- gsvaParam(exprData = expr, geneSets = modules)
  gsva.modules <- as.data.frame(gsva(gsvaPar, verbose=FALSE))
  
  return(gsva.modules)
}

calc_spongeffects_ssgsea <- function(expr, modules, n_cores=1) {
  
  ssgseaPar <- ssgseaParam(exprData = expr, geneSets = modules)
  ssgsea.modules <- as.data.frame(gsva(ssgseaPar, verbose=FALSE))
  
  return(ssgsea.modules)
}


create_Sponge_modules <- function(var, ceRNA_interactions, ceRNA_centralities, mscor=0.1, padj=0.05,
                                  feature_col='feature_type', RNAs=c("lncRNA","protein_coding"), ensembl_col="ensembl_id", max_modules, n_cores) {
  
  if(! feature_col %in% colnames(var)) {
    stop("You chose an invalid column for the type of RNA.")
  }
  
  if(ensembl_col %in% colnames(var)) {
    RNAs.ofInterest <- var[get(feature_col) %in% c("lncRNA", "protein_coding"), get(ensembl_col)]
  } else {
    stop("You chose an invalid column for the ensembl ids.")
  }
  
  # filter network
  filtered_network_centralities <- filter_ceRNA_network(sponge_effects = ceRNA_interactions, Node_Centrality = ceRNA_centralities, 
                                                        add_weighted_centrality=TRUE, mscor.threshold = mscor, padj.threshold = padj)
  
  central_gene_modules <- get_central_modules(central_nodes = RNAs.ofInterest, node_centrality = filtered_network_centralities$Node_Centrality, 
                                            ceRNA_class = RNAs, centrality_measure = "Weighted_Degree", cutoff = max_modules)
  
  Sponge.modules <- define_modules(network = filtered_network_centralities$Sponge.filtered, central.modules = central_gene_modules, remove.central = FALSE, set.parallel = n_core > 1)
  
  # filter modules to at least contain 10 ids
  Sponge.modules <- Sponge.modules[sapply(Sponge.modules, length) >= 10]
  
  return(Sponge.modules)
}


# for AUCell and SPONGeffects
create_Genie_modules <- function(regulon_df, expr, top=100000, dir, k, n_regulons) {
  # filter Top n percent of edges
  
  regulon_df %>%
    arrange(weight) %>%
    head(top)
  
  # calculate regulons
  setnames(regulon_df, c("regulatoryGene", "targetGene", "weight"), c("tf", "target", "mi"))
  
  file <- file.path(dir, "genie3_aracne_format.txt")
  fwrite(regulon_df, file, sep = "\t", col.names = FALSE)
  
  setnames(regulon_df, c("tf", "target", "mi"), c("regulatoryGene", "targetGene", "weight"))
  
  expr_filt <- expr[apply(expr, 1, sd) > 0, ]
  regulons <- aracne2regulon(afile = file, eset = expr_filt, format = "3col", verbose = FALSE)
  
  # filter regulons to at least contain 10 ids
  regulons <- regulons[sapply(regulons, function(r) length(r$tfmode) >= 10)]
  
  # calculate top n_regulons by mean weight of top k target genes
  regulons <- regulons[order(sapply(regulons, mean_top_n_likelihood, k), decreasing = TRUE)]
  regulons <- head(regulons, n_regulons)
  
  return(regulons)
}


mean_top_n_likelihood <- function(regulon, n) {
  weights <- regulon$likelihood
  top_weights <- sort(weights, decreasing = TRUE)[1:min(n, length(weights))]
  return(mean(top_weights))
}



filter_network <- function(network, genesets, type=c("Sponge", "Genie")) {
  genes <- unique(c(unlist(genesets), names(genesets)))
  
  if (type == "Sponge") {
    net_filt <- network[(geneA %in% genes) | (geneB %in% genes)]
  } else if (type == "Genie") {
    net_filt <- network[(regulatoryGene %in% genes) | (targetGene %in% genes)]
  } else {
    net_filt <- NA
  }
  
  return(net_filt)
}
