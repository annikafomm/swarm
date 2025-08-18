library(SPONGE)
library(dplyr)
library(data.table)
library(tibble)
library(GSVA)


# file: "../datasets_prepro_leidenalg/GSM6592055_M8-prepro.h5ad"
# "../networks/SPONGE/breast_invasive_carcinoma/breast_invasive_carcinoma_interactionNetwork.csv"
# "../networks/SPONGE/breast_invasive_carcinoma/breast_invasive_carcinoma_networkAnalysis.csv"


# TODO: choose parameters wisely
calc_spongeEffects <- function(expr, modules, oe=TRUE, gsva=TRUE, ssgsea=TRUE, n_cores=1) {
  
  modules <- modules[sapply(modules, length) > 1]
  
  oe.modules <- NA
  if (oe) {
    oe.modules <- enrichment_modules(Expr.matrix = expr, modules = modules, bin.size = 10, min.size = 1, max.size = Inf, min.expr = 1, method = "OE", cores=n_cores)
  }
  
  gsva.modules <- NA
  if (gsva) {
    gsvaPar <- gsvaParam(exprData = expr, geneSets = modules)
    gsva.modules <- as.data.frame(gsva(gsvaPar, verbose=FALSE))
  }
  
  ssgsea.modules <- NA
  if (ssgsea) {
    ssgseaPar <- ssgseaParam(exprData = expr, geneSets = modules)
    ssgsea.modules <- as.data.frame(gsva(ssgseaPar, verbose=FALSE))
  }
  
  return(list(oe = oe.modules,
              gsva = gsva.modules,
              ssgsea = ssgsea.modules))
}


create_Sponge_modules <- function(var, ceRNA_interactions, ceRNA_centralities, mscor=0.1, padj=0.05,
                                  feature_col='feature_type', RNAs=c("lncRNA","protein_coding"), ensembl_col="V1", n_modules, n_cores) {
  
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
                                            ceRNA_class = RNAs, centrality_measure = "Weighted_Degree", cutoff = n_modules)
  
  Sponge.modules <- define_modules(network = filtered_network_centralities$Sponge.filtered, central.modules = central_gene_modules, remove.central = FALSE, set.parallel = n_core > 1)
  
  return(Sponge.modules)
}


# for AUCell and SPONGeffects
create_Genie_modules <- function(regulon_df, weight.threshold, k, cutoff) {
  # list of regulons (= named list of weights, names are regulators)
  genie_network <- filter(regulon_df, weight > weight.threshold)
  
  # calculate top 'cutoff' regulons by mean weight of top k target genes
  top_regulons <- genie_network %>%
    group_by(regulatoryGene) %>%
    slice_max(order_by = weight, n = k) %>%
    summarize(mean_weight = mean(weight, na.rm = TRUE)) %>%
    arrange(desc(mean_weight), by_group=TRUE) %>%
    head(cutoff) %>% 
    pull(regulatoryGene)
  
  regulons <- genie_network %>%
    filter(regulatoryGene %in% top_regulons) %>%
    { split(.$targetGene, .$regulatoryGene) }
  
  return(regulons)
}
