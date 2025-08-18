library(Matrix)
library(data.table)

source("./calc_spongEffects.R")

library(optparse)

main <- function() {
  option_list <- list(
    # output dir
    make_option("--outdir", type = "character", help = "Output directory"),
    
    # network files
    make_option("--sponge_network", type="character", help="Path to Sponge interactionNetwork file", metavar="file"),
    make_option("--sponge_analysis", type="character", help="Path to Sponge networkAnalysis file", metavar="file"),
    make_option("--genie_network", type="character", help="Path to GENIE3 network file", metavar="file"),
    
    # expression files
    make_option("--expr_matrix", type="character", help="Path to expression matrix file", metavar="file"),
    make_option("--var", type="character", help="Path to var file", metavar="file"),
    make_option("--cells", type="character", help="Path to cells file", metavar="file"),
    make_option("--genes", type="character", help="Path to genes file", metavar="file"),
    
    # score flags
    make_option("--aucell", action="store_true", default=FALSE, help="Calculate AUCell score"),
    make_option("--oe", action="store_true", default=FALSE, help="Calculate SPONGeffect score: OE"),
    make_option("--gsva", action="store_true", default=FALSE, help="Calculate SPONGeffect score: GSVA"),
    make_option("--ssgsea", action="store_true", default=FALSE, help="Calculate SPONGeffect score: ssGSEA"),
    make_option("--viper", action="store_true", default=FALSE, help="Calculate VIPER score"),
    
    # Sponge modules
    make_option("--mscor", type = "numeric", default = 0.1, help = "mscor threshold", metavar = "number"),
    make_option("--padj", type = "numeric", default = 0.05, help = "padj threshold", metavar = "number"),
    make_option("--feature_col", type = "character", default = 'feature_type', help = "Column name with type of RNA", metavar = "character"),
    make_option("--ensembl_col", type = "character", default = 'V1', help = "Column name with ensembl ids", metavar = "character"),
    make_option("--RNA_types", type = "character", default = "lncRNA,protein_coding", help = "Comma-separated list of feature types", metavar = "types"),
    make_option("--n_modules", type = "numeric", default = 20, help = "maximum number of modules", metavar = "number"),
    make_option("--n_cores", type = "numeric", default = 1, help = "number of cores", metavar = "number"),
    
  )
  
  parser <- OptionParser(option_list=option_list)
  args <- parse_args(parser)
  
  # output dir
  
  if (is.null(args.outdir)) {
    stop("Please provide an output dir.")
  }
  if (!dir.exists(args$outdir)) dir.create(args$outdir, recursive = TRUE)
  
  
  # required flags
  
  if ((is.null(args$sponge_network) | is.null(args$sponge_analysis)) & is.null(args$genie_network)) {
    stop("Please provide a SPONGE or GENIE3 network with --sponge_network and --sponge_analysis or --genie_network.")
  }
  if (is.null(args$expr_matrix)) {
    stop("Please provide an expression matrix with --expr_matrix.")
  }
  
  if (!(is.null(args$sponge_network) | is.null(args$sponge_analysis)) & is.null(args$var)) {
    stop("Please provide a var df with --var.")
  }
  
  if (is.null(args$var) & is.null(args$genes)) {
    stop("Please provide a genes list with --genes or a var df with --var.")
  }
  
  if (! (args.aucell | args.oe | args.gsva | args.ssgsea | args.viper)) {
    stop("Please choose a score to calculate.")
  }
  
  
  # expression matrix
  
  expr <- readMM(args.expr_matrix)
  if (! is.null(args.cells)) {
    rownames(expr) <- fread(args.cells, header = FALSE)$V1
  }
  if (is.null(args.var)) {
    colnames(expr) <- fread(args.genes, header = FALSE)$V1
  } else {
    var_df <- fread(args.var)
    colnames(expr) <- var_df[,get(args.ensembl_col)]
  }
  expr <- t(as.matrix(expr))
  
  
  # sponge modules
  
  if (! (is.null(args$sponge_network) | is.null(args$sponge_analysis)) & (args.aucell | args.oe | args.gsva | args.ssgsea)) {
    ceRNA_interactions <- fread(args.sponge_network)
    ceRNA_centralities <- fread(args.sponge_analysis)
    
    sponge_modules <- create_Sponge_modules(var_df, ceRNA_interactions, ceRNA_centralities, mscor=args.mscor, padj=args.padj,
                                      feature_col=args.feature_col, RNAs=args.RNA_types, ensembl_col=args.ensembl_col, n_modules=args.n_modules, n_cores=args.n_cores)
    
    
    if (args.oe | args.gsva | args.ssgsea) {
      sponge_effects <- calc_spongeEffects(expr, sponge_modules, oe=args.oe, gsva=args.gsva, ssgsea=args.ssgsea, n_cores=args.n_cores)
      
      
      file_clean <- sub("_X$", "", tools::file_path_sans_ext(basename(args.expr_matrix)))
      lapply(names(effects), function(x) file.path(args.outdir, paste0("spongeffects_", x, "_sponge_", file_clean, ".csv")))
    }
    
    if (args.aucell) {
      # TODO
    }
    
    
  }
  
  
  
  # TODO: regulons for sponge/aucell and viper
  
  if (! is.null(args$genie_network)) {
    regulon_df <- fread(args.genie_network)
  }
  
  
  
  
  
}


# only run main() if the script is executed directly, not sourced
if (identical(environment(), globalenv()) && !length(sys.frames())) {
  main()
}
