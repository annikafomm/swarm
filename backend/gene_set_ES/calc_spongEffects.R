library(SPONGE)
library(doParallel)
library(foreach)
library(dplyr)
library(anndata)
library(data.table)
library(tibble)

# Register your backend here
num.of.cores <- 1
cl <- makeCluster(num.of.cores) 
registerDoParallel(cl)

# file: "../datasets_prepro_leidenalg/GSM6592055_M8-prepro.h5ad"
# "../networks/SPONGE/breast_invasive_carcinoma/breast_invasive_carcinoma_interactionNetwork.csv"
# "../networks/SPONGE/breast_invasive_carcinoma/breast_invasive_carcinoma_networkAnalysis.csv"


# TODO: choose parameters wisely
calc_spongeEffects <- function(adata, interaction_net, network_analysis, 
                               mscor=0.1, padj=0.05, 
                               feature_col='feature_type', RNAs=c("lncRNA","protein_coding"), 
                               ensembl_col=NA, n_modules=20, 
                               oe=TRUE, gsva=TRUE, ssgsea=TRUE, n_cores=1) {
  
  # Expression matrix (may be sparse)
  expr <- adata$X
  
  # Convert to dense matrix if needed
  if (inherits(expr, "dgRMatrix")) {
    expr <- as.matrix(expr)
  }
  
  expr <- t(expr)
  
  
  # TODO: How to do that universal?
  if (is.na(ensembl_col)) {
    RNAs.ofInterest <- adata$var %>% 
      filter(feature_col %in% RNAs) %>%
      rownames()
  } else if(ensembl_col %in% adata$var_keys()) {
    RNAs.ofInterest <- adata$var %>% 
      filter(feature_col %in% RNAs) %>%
      pull(ensembl_col)
  } else {
    print("You chose an invalid column for the ensembl ids.")
  }
  print(length(RNAs.ofInterest))
  
  create_Sponge_modules <- function(interaction_net, network_analysis, mscor, padj,
                                    RNAs, RNAs.ofInterest, n_modules, n_cores)
  
  oe.modules <- None
  if (oe) {
    oe.modules <- enrichment_modules(Expr.matrix = expr, modules = Sponge.modules, bin.size = 10, min.size = 1, max.size = Inf, min.expr = 1, method = "OE", cores=n_cores)
  }
  
  gsva.modules <- None
  if (gsva) {
    gsvaPar <- gsvaParam(exprData = expr, geneSets = Sponge.modules)
    gsva.modules <- as.data.frame(gsva(gsvaPar, verbose=FALSE))
  }
  
  ssgsea.modules <- None
  if (ssgsea) {
    ssgseaPar <- ssgseaParam(exprData = expr, geneSets = Sponge.modules)
    ssgsea.modules <- as.data.frame(gsva(ssgseaPar, verbose=FALSE))
  }
  
  return (oe.modules, gsva.modules, ssgsea.modules)
}


create_Sponge_modules <- function(interaction_net, network_analysis, mscor, padj,
                                  RNAs, RNAs.ofInterest, n_modules, n_cores) {
  ceRNA_interactions <- fread(interaction_net)
  ceRNA_centralities <- fread(network_analysis)
  
  # filter network
  filtered_network_centralities <- filter_ceRNA_network(sponge_effects = ceRNA_interactions, Node_Centrality = ceRNA_centralities, 
                                                        add_weighted_centrality=T, mscor.threshold = mscor, padj.threshold = padj)
  
  central_gene_modules<-get_central_modules(central_nodes = RNAs.ofInterest, node_centrality = filtered_network_centralities$Node_Centrality, 
                                            ceRNA_class = RNAs, centrality_measure = "Weighted_Degree", cutoff = n_modules)
  
  
  
  Sponge.modules <- define_modules(network = filtered_network_centralities$Sponge.filtered, central.modules = central_gene_modules, remove.central = F, set.parallel = n_core > 1)
  
  return(Sponge.modules)
}


create_Genie_modules <- function() {
  # Load the regulon table
  regulon_df <- fread("../networks/GENIE3/BRCA/genie3_BRCA_mrn.top_100k.csv")
  
  # list of regulons (= named list of weights, names are regulators)
  regulon_list <- regulon_df %>%
    group_by(regulatoryGene) %>%
    summarize(regulon = list(setNames(weight, targetGene))) %>%
    deframe()
}