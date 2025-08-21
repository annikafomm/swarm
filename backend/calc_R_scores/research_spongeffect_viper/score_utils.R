library(data.table)
library(tibble)


filter_genie_network <- function(genie_network, weight.threshold) {
   filter(genie_network, weight > weight.threshold)
}


get_top_regulons <- function(genie_network, k, cutoff = 1000) {
  top_regulons <- genie_network %>%
    group_by(regulatoryGene) %>%
    slice_max(order_by = weight, n = k) %>%
    summarize(mean_weight = mean(weight, na.rm = TRUE)) %>%
    arrange(desc(mean_weight), by_group=TRUE) %>%
    head(cutoff) %>%
    pull(regulatoryGene)
  
  genie_network %>%
    filter(regulatoryGene %in% top_regulons) %>%
    { split(.$targetGene, .$regulatoryGene) }
}
















