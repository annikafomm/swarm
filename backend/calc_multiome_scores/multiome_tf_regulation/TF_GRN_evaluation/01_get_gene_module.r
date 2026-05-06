
# GET GENE MODULE FOR GO TERM
# conda activate mzb_R
# Rscript 01_get_gene_module.r GO:0030199 tmp
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
    stop("Usage: Rscript script.r <GO_id> <outdir>")
}

GO_id <- args[1]
outdir <- args[2]

library(AnnotationDbi)
library(org.Hs.eg.db)
gene_module <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = GO_id,
  keytype = "GOALL",
  columns = c("ENSEMBL", "SYMBOL", "GOALL", "ONTOLOGYALL")
)

# write to file for use in python script
write.csv(gene_module, file = paste0(outdir, "/gene_module.csv"), row.names = FALSE)
