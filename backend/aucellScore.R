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
                                            ceRNA_class = RNAs, centrality_measure ="Weighted_Degree", cutoff = n_modules)

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

convert_ensembl_to_symbols <- function(df, col) {
  if (missing(col) || !col %in% colnames(df)) {
    stop(paste("Column", col, "not found in the dataframe."))
  }
  
  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    stop("Package 'biomaRt' is required but not installed.")
  }
  
  # Connect to Ensembl mart
  mart <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  
  # Get gene symbols for Ensembl IDs
  gene_symbols <- biomaRt::getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol"),
    filters = "ensembl_gene_id",
    values = df[[col]],
    mart = mart
  )
  
  # Merge gene symbols into original df (inner join to keep only matches)
  df_merged <- merge(df, gene_symbols, by.x = col, by.y = "ensembl_gene_id", all.x = FALSE, sort = FALSE)
  
  # Remove rows with missing or empty gene symbols
  df_merged <- df_merged[!(is.na(df_merged$hgnc_symbol) | df_merged$hgnc_symbol == ""), ]
  
  # Replace Ensembl IDs with gene symbols in the original column name
  df_merged[[col]] <- df_merged$hgnc_symbol
  
  # Drop the temporary hgnc_symbol column
  df_merged$hgnc_symbol <- NULL
  
  # Return the updated dataframe
  return(df_merged)
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

ad <- zellkonverter::readH5AD(adata)
interaction_network <- fread(interaction_network, strip.white=TRUE)
network_analysis <- fread(network_analysis, strip.white=TRUE)

if (!is.null(interaction_network) && !is.null(network_analysis)) {
  RNAs.ofInterest <- calc_spongeEffects(ad, mscor = 0.3, padj = 0.05, 
                                        feature_col = 'feature_type', 
                                        RNAs = c("lncRNA", "protein_coding"), 
                                        ensembl_col = ensembl_col, 
                                        n_modules = 20, 
                                        oe = TRUE, gsva = TRUE, ssgsea = TRUE, 
                                        n_cores = parallel::detectCores())

print(paste("Identified", length(RNAs.ofInterest), "genes of interest from the AnnData object."))


#interaction_network <- convert_ensembl_to_symbols(interaction_network, 'geneA')
#interaction_network <- convert_ensembl_to_symbols(interaction_network, 'geneB')
#network_analysis <- convert_ensembl_to_symbols(network_analysis, 'geneB')

# Load the SPONGE network
sponge_network <- create_Sponge_modules(interaction_network, network_analysis, mscor = 0.03, padj = 0.05,
                                        RNAs = c("ceRNA", "lncRNA", "mRNA"),
                                        RNAs.ofInterest = RNAs.ofInterest,
                                        n_modules = 10, n_cores = parallel::detectCores())

# Calculate AUCell scores
aucell_scores_sponge <- aucell_score(sponge_network, ad, ensembl_col)

}

if (!is.null(genie_3_network)) {
  # Create gene sets from the Genie3 network
  genie3_gene_sets <- create_Genie3_gene_sets(genie_3_network, n_cores = parallel::detectCores())

  aucell_scores_genie3 <- aucell_score(genie3_gene_sets, ad, "")
}

# Get dataframes from AUCell scores

aucell_scores_sponge <- AUCell::getAUC(aucell_scores_sponge)
aucell_scores_sponge <- as.data.frame(aucell_scores_sponge)

aucell_scores_genie3 <- AUCell::getAUC(aucell_scores_genie3)
aucell_scores_genie3 <- as.data.frame(aucell_scores_genie3)

adata_prefix <- sub(".h5ad", "", basename(adata))

  
# Save to csv
write.csv(aucell_scores_sponge, file = paste0("aucell_scores_sponge_", adata_prefix, ".csv"), row.names = TRUE)
write.csv(aucell_scores_genie3, file = paste0("aucell_scores_genie3_", adata_prefix, ".csv"), row.names = TRUE)





