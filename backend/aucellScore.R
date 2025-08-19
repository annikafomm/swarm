library('SPONGE')
library('zellkonverter')
library('SingleCellExperiment')
library('dplyr')
library('data.table')
library('tibble')
library('biomaRt')
library('AUCell')


calc_spongeEffects <- function(adata,
                               mscor=0.1, padj=0.05,
                               feature_col='feature_type', RNAs=c("lncRNA","protein_coding"),
                               ensembl_col=NA, n_modules=20,
                               oe=TRUE, gsva=TRUE, ssgsea=TRUE, n_cores=1) {

  # Get Expresision matrix from SingleCellExperiment object
  expr <- assay(adata, "X")
  expr <- as.matrix(expr)
  expr <- t(expr)
  # Get var (genes)
  var <- rowData(adata)

  if (is.na(ensembl_col)) {
      RNAs.ofInterest <- rownames(var[as.character(var[[feature_col]]) %in% RNAs, ])
  } else if (ensembl_col %in% colnames(var)) {
      RNAs.ofInterest <- var[as.character(var[[feature_col]]) %in% RNAs, ensembl_col]
  } else {
      message("You chose an invalid column for the Ensembl IDs.")
      RNAs.ofInterest <- character(0)
  }

  return(RNAs.ofInterest)

}

create_Sponge_modules <- function(interaction_net, network_analysis, mscor, padj,
                                  RNAs, RNAs.ofInterest, n_modules, n_cores) {
  # Filter by mscor and padj
  filtered_network_centralities <- filter_ceRNA_network(sponge_effects = interaction_net, Node_Centrality = network_analysis,
                                                        add_weighted_centrality=TRUE, mscor.threshold = mscor, padj.threshold = padj)
  print(filtered_network_centralities$Node_Centrality)

  central_gene_modules<-get_central_modules(central_nodes = RNAs.ofInterest, node_centrality = filtered_network_centralities$Node_Centrality,
                                            ceRNA_class = RNAs, centrality_measure ="Weighted_Degree", cutoff = 5)

  intersect_genes <- intersect(RNAs.ofInterest, filtered_network_centralities$Node_Centrality$gene)

  print(paste("Number of central genes:", length(central_gene_modules)))
  print(paste("Number of genes in the network:", length(intersect_genes)))


  Sponge.modules <- define_modules(network = filtered_network_centralities$Sponge.filtered, central.modules = central_gene_modules, remove.central = FALSE, set.parallel = n_cores > 1)

  return(Sponge.modules)
}

create_Genie3_gene_sets <- function(genie3network, n_cores, k = 100, cutoff = 50) {
  # Load the regulon table
  regulon_df <- fread(genie3network)

  regulon_df <- filter(regulon_df, weight > 0.01)

  top_regulons <- regulon_df %>%
    group_by(regulatoryGene) %>%
    slice_max(order_by = weight, n = k) %>%
    summarize(mean_weight = mean(weight, na.rm = TRUE)) %>%
    arrange(desc(mean_weight), by_group=TRUE) %>%
    head(cutoff) %>%
    pull(regulatoryGene)

  # Create gene sets for AUCell
  regulons_for_aucell <- regulon_df %>%
    filter(regulatoryGene %in% top_regulons) %>%
    group_by(regulatoryGene) %>%
    summarize(targets = list(targetGene), .groups = 'drop') %>%
    deframe()

  return(regulons_for_aucell)

}



aucell_score <- function(sponge_network, anndata, ensembl_col) {

  # Extract gene expression matrix
  expr_matrix <- assay(anndata, "X")



  if (!is.null(ensembl_col) && ensembl_col != "") {
    # Get the Ensembl IDs from the AnnData object
    if (!ensembl_col %in% colnames(rowData(anndata))) {
      stop("Ensembl column not found in the AnnData object.")
    }
    ensembl_ids <- rowData(anndata)[[ensembl_col]]

  # Ensure the row names of the expression matrix are Ensembl IDs
    if (length(ensembl_ids) != nrow(expr_matrix)) {
      stop("The number of Ensembl IDs does not match the number of rows in the expression matrix.")
    }
    rownames(expr_matrix) <- ensembl_ids
  }

  # Create AUCell object
  rankings <- AUCell::AUCell_buildRankings(expr_matrix)

  # Calculate AUCell scores for the SPONGE network
  aucell_scores <- AUCell::AUCell_calcAUC(sponge_network, rankings, aucMaxRank = nrow(expr_matrix), nCores = parallel::detectCores())

  return(aucell_scores)

}


args = commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: Rscript backend/aucellScore.R <interaction_network> <network_analysis> <genie_3_network> <anndata> <ensembl_col>")

}


interaction_network <- args[1]
network_analysis <- args[2]
genie_3_network <-  args[3]
adata <-args[4]
ensembl_col <-  args[5]

adata_prefix <- sub(".h5ad", "", basename(adata))

ad <- zellkonverter::readH5AD(adata)

interaction_network <- fread(interaction_network, strip.white=TRUE)
network_analysis <- fread(network_analysis, strip.white=TRUE)

if (!is.null(interaction_network) && !is.null(network_analysis)) {
  RNAs.ofInterest <- calc_spongeEffects(ad, mscor = 0.05, padj = 0.01,
                                        feature_col = 'feature_type',
                                        RNAs = c("lncRNA", "protein_coding"),
                                        ensembl_col = ensembl_col,
                                        n_modules = 10,
                                        oe = TRUE, gsva = TRUE, ssgsea = TRUE,
                                        n_cores = parallel::detectCores())

print(paste("Identified", length(RNAs.ofInterest), "genes of interest from the AnnData object."))

# Load the SPONGE network
sponge_network <- create_Sponge_modules(interaction_network, network_analysis, mscor = 0.03, padj = 0.05,
                                        RNAs = c("ceRNA", "lncRNA", "mRNA"),
                                        RNAs.ofInterest = RNAs.ofInterest,
                                        n_modules = 10, n_cores = parallel::detectCores())


# Filter Sponge network by genes that are in the AnnData object
sponge_network <- sponge_network[names(sponge_network) %in% RNAs.ofInterest]
sponge_network <- lapply(sponge_network, function(gene_set) gene_set[gene_set %in% RNAs.ofInterest])
# Save to json file
sponge_genesets_file <- paste0("sponge_gene_sets_", adata_prefix, ".json")
jsonlite::write_json(sponge_network, path = sponge_genesets_file, pretty = TRUE)

print(interaction_network)
# Filter interaction network by adjusted p-value and mscor
filtered_network <- interaction_network[interaction_network$p.adj < 0.05 & abs(interaction_network$mscor) > 0.03, ]
filtered_network <- filtered_network[filtered_network$geneA %in% RNAs.ofInterest & filtered_network$geneB %in% RNAs.ofInterest, ]

# Save the filtered network to a file
filtered_network_file <- paste0("breast_invasive_carcinoma_filtered_Network.tsv")
write.table(filtered_network, file = filtered_network_file, sep = "\t", quote = FALSE, row.names = FALSE)

# Calculate AUCell scores
aucell_scores_sponge <- aucell_score(sponge_network, ad, ensembl_col)
aucell_scores_sponge <- AUCell::getAUC(aucell_scores_sponge)
aucell_scores_sponge <- as.data.frame(aucell_scores_sponge)
write.csv(aucell_scores_sponge, file = paste0("aucell_scores_sponge_", adata_prefix, ".csv"), row.names = TRUE)
}

if (!is.null(genie_3_network)) {
  # Create gene sets from the Genie3 network
  genie_3_gene_sets_file <- paste0("genie3_network_", adata_prefix, ".txt")
  genie3_gene_sets <- create_Genie3_gene_sets(genie_3_network, n_cores = parallel::detectCores())

  # Save the Genie3 as a JSON file
  genie3_gene_sets_file <- paste0("genie3_gene_sets_", adata_prefix, ".json")
  jsonlite::write_json(genie3_gene_sets, path = genie3_gene_sets_file, pretty = TRUE)

  aucell_scores_genie3 <- aucell_score(genie3_gene_sets, ad, NULL)
  aucell_scores_genie3 <- AUCell::getAUC(aucell_scores_genie3)
  aucell_scores_genie3 <- as.data.frame(aucell_scores_genie3)
  write.csv(aucell_scores_genie3, file = paste0("aucell_scores_genie3_", adata_prefix, ".csv"), row.names = TRUE)
}
