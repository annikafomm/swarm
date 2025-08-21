suppressPackageStartupMessages({
  suppressWarnings(library(Matrix))
  suppressWarnings(library(data.table))
  suppressWarnings(library(viper))
  suppressWarnings(library(jsonlite))
  suppressWarnings(library(optparse))
  suppressWarnings(library(SPONGE))
  suppressWarnings(library(dplyr))
  suppressWarnings(library(tibble))
  suppressWarnings(library(GSVA))
  suppressWarnings(library(AUCell))
})

source("./utils.R")

#  Rscript calc_scores.R   --dir ../datasets_scores/GSM6592049_M2_prepro   --dir ../datasets_scores/GSM6592049_M2_scores   --sponge_network ../networks/SPONGE/breast_invasive_carcinoma/breast_invasive_carcinoma_interactionNetwork.csv   --sponge_analysis ../networks/SPONGE/breast_invasive_carcinoma/breast_invasive_carcinoma_networkAnalysis.csv --genie_network ../networks/GENIE3/BRCA/genie3_BRCA_mrn.top_100k.csv --ensembl_col ensemble_id  --aucell --oe --gsva --ssgsea --viper

# ! parallel::detectCores()

# logging helper
log_message <- function(msg, logfile, indent = 0) {
  prefix <- paste(rep(" ", indent), collapse = "")
  line <- paste0(prefix, msg)   # <-- no "\n" here
  cat(line, "\n", sep = "")     # console
  write(line, file = logfile, append = TRUE)  # file (no extra newline)
}

format_runtime <- function(t0) {
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  mins <- floor(elapsed / 60)
  secs <- round(elapsed %% 60)
  sprintf("%d min %d sec", mins, secs)
}

main <- function() {
  option_list <- list(
    # dir with files
    make_option("--dir", type="character", help="Dir containing all input files (expr.mtx, var.tsv, (cells.txt))"),
    make_option("--log", type="character", help="Path to the log file", metavar="file"),
    
    # network files
    make_option("--sponge_network", type="character", help="Path to Sponge interactionNetwork file", metavar="file"),
    make_option("--sponge_analysis", type="character", help="Path to Sponge networkAnalysis file", metavar="file"),
    make_option("--genie_network", type="character", help="Path to GENIE3 network file", metavar="file"),
    
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
    make_option("--ensembl_col", type = "character", default = 'ensembl_id', help = "Column name with ensembl ids", metavar = "character"),
    make_option("--RNA_types", type = "character", default = "lncRNA,protein_coding", help = "Comma-separated list of feature types", metavar = "types"),
    make_option("--max_modules", type = "numeric", default = 20, help = "maximum number of modules", metavar = "number"),
    make_option("--n_cores", type = "numeric", default = 1, help = "number of cores", metavar = "number"),
    
    # regulons
    # TODO: best practice cutoff for weigths and k_reg_genes
    make_option("--top_n", type = "numeric", default = 100000, help = "keep top n weights of genie network", metavar = "number"),
    make_option("--k_reg_genes", type = "numeric", default = 20, help = "number of regulatory genes for calculating the mean weight for identifying the top n regulons", metavar = "number"),
    make_option("--n_regulons", type = "numeric", default = 20, help = "number of regulons", metavar = "number")
    
  )
  
  parser <- OptionParser(option_list=option_list)
  args <- parse_args(parser)
  
  if (is.null(args$dir)) {
    stop("Please provide an output dir.")
  }
  if (!dir.exists(args$dir)) dir.create(args$dir, recursive = TRUE)
  
  if (! is.null(args$log)) logfile <- args$log
  log_message(paste0("Pipeline started at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n"), logfile)
  
  log_message("Input parameters (JSON):", logfile)
  log_message(toJSON(args, pretty = TRUE, auto_unbox = TRUE), logfile)
  
  # directories
  if (is.null(args$dir)) {
    stop("Please provide an input dir.")
  }
  if (! all(c("expr.mtx", "var.tsv") %in% list.files(args$dir))) {
    stop("The input dir you provided does not contain the right files (expr.mtx & var.tsv)")
  }
  
  # required flags
  
  if ((is.null(args$sponge_network) | is.null(args$sponge_analysis)) & is.null(args$genie_network)) {
    stop("Please provide a SPONGE or GENIE3 network with --sponge_network and --sponge_analysis or --genie_network.")
  }
  
  if (! (args$aucell | args$oe | args$gsva | args$ssgsea | args$viper)) {
    stop("Please choose a score to calculate.")
  }
  
  
  # expression matrix
  log_message("Loading expression matrix (and additional files) ...", logfile)
  t0 <- Sys.time()
  
  expr <- readMM(file.path(args$dir, "expr.mtx"))
  if (file.exists(file.path(args$dir, "cells.txt"))) {
    rownames(expr) <- fread(file.path(args$dir, "cells.txt"), header = FALSE)$V1
  }
  var_df <- fread(file.path(args$dir, "var.tsv"))
  expr <- t(as.matrix(expr))
  
  if (nrow(expr) != nrow(var_df)) {
    stop("Expression matrix is not as expected.")
  }
  
  log_message(sprintf("Expression matrix loaded in %s", format_runtime(t0)), logfile, 2)
  
  
  
  # SPONGE networks
  
  if (! is.null(args$sponge_network) & !is.null(args$sponge_analysis) & (args$aucell | args$oe | args$gsva | args$ssgsea)) {
    log_message("Loading SPONGE network and create SPONGE modules ...", logfile)
    t0 <- Sys.time()
    
    # load files and calculate SPONGE modules
    ceRNA_interactions <- fread(args$sponge_network)
    ceRNA_centralities <- fread(args$sponge_analysis)
    
    log_message(sprintf("SPONGE network loaded in %s", format_runtime(t0)), logfile, 2)
    
    rownames(expr) <- var_df[,get(args$ensembl_col)]
    
    t0 <- Sys.time()
    
    sponge_modules <- create_Sponge_modules(var_df, ceRNA_interactions, ceRNA_centralities, mscor=args$mscor, padj=args$padj,
                                      feature_col=args$feature_col, RNAs=args$RNA_types, ensembl_col=args$ensembl_col, max_modules=args$max_modules, n_cores=args$n_cores)
    
    log_message(sprintf("SPONGE modules created in %s", format_runtime(t0)), logfile, 2)
    
    modules_exist <- TRUE
    if (length(names(sponge_modules)) == 0) {
      log_message("With the chosen parameters no SPONGE modules could be created.", logfile, 2)
      modules_exist <- FALSE
    } else {
      t0 <- Sys.time()
      
      write_json(sponge_modules, path = file.path(args$dir, paste0("sponge_genesets.json")), pretty = TRUE, auto_unbox = TRUE)
      
      log_message(sprintf("SPONGE modules written in %s", format_runtime(t0)), logfile, 2)
    }
    
    # calculate SPONGeffects scores
    if (args$oe & modules_exist) {
      log_message("Calculating SPONGeffects-OE scores (SPONGE) ...", logfile)
      t0 <- Sys.time()
      
      oe_scores <- calc_spongeffects_oe(expr, sponge_modules, n_cores=args$n_cores)
      
      log_message(sprintf("SPONGeffects-OE scores computed in %s", format_runtime(t0)), logfile, 2)
      t0 <- Sys.time()
      
      write.csv(oe_scores, file = file.path(args$dir, paste0("spongeffects_OE_scores_sponge.csv")), row.names = TRUE)
      
      log_message(sprintf("SPONGeffects-OE scores written in %s", format_runtime(t0)), logfile, 2)
    }
    
    if (args$gsva & modules_exist) {
      log_message("Calculating SPONGeffects-GSVA scores (SPONGE) ...", logfile)
      t0 <- Sys.time()
      
      gsva_scores <- calc_spongeffects_gsva(expr, sponge_modules, n_cores=args$n_cores)
      
      log_message(sprintf("SPONGeffects-GSVA scores computed in %s", format_runtime(t0)), logfile, 2)
      t0 <- Sys.time()
      
      write.csv(gsva_scores, file = file.path(args$dir, paste0("spongeffects_GSVA_scores_sponge.csv")), row.names = TRUE)
      
      log_message(sprintf("SPONGeffects-GSVA scores written in %s", format_runtime(t0)), logfile, 2)
    }
    
    if (args$ssgsea & modules_exist) {
      log_message("Calculating SPONGeffects-ssGSEA scores (SPONGE) ...", logfile)
      t0 <- Sys.time()
      
      ssgsea_scores <- calc_spongeffects_ssgsea(expr, sponge_modules, n_cores=args$n_cores)
      
      log_message(sprintf("SPONGeffects-ssGSEA scores computed in %s", format_runtime(t0)), logfile, 2)
      t0 <- Sys.time()
      
      write.csv(ssgsea_scores, file = file.path(args$dir, paste0("spongeffects_ssGSEA_scores_sponge.csv")), row.names = TRUE)
      
      log_message(sprintf("SPONGeffects-ssGSEA scores written in %s", format_runtime(t0)), logfile, 2)
    }
    
    # calculate AUCell scores
    if (args$aucell & modules_exist) {
      log_message("Calculating AUCell scores (SPONGE) ...", logfile)
      t0 <- Sys.time()
      
      aucell_scores <- calc_aucell_score(expr, sponge_modules, args$n_cores)
      
      log_message(sprintf("AUCell scores computed in %s", format_runtime(t0)), logfile, 2)
      t0 <- Sys.time()
      
      write.csv(aucell_scores, file = file.path(args$dir, paste0("aucell_scores_sponge.csv")), row.names = TRUE)
      
      log_message(sprintf("AUCell scores written in %s", format_runtime(t0)), logfile, 2)
    }
  }
  
  
  # GENIE3 networks
  # TODO: regulons for sponge/aucell and viper
  
  if (! is.null(args$genie_network)) {
    log_message("Loading GENIE3 network and create regulons ...", logfile)
    
    rownames(expr) <- var_df[,V1]
    
    t0 <- Sys.time()
    # load network file and calculate regulons
    regulon_df <- fread(args$genie_network)
    
    log_message(sprintf("GENIE3 network loaded in %s", format_runtime(t0)), logfile, 2)
    t0 <- Sys.time()
    
    vregulons <- create_Genie_modules(regulon_df, expr=expr, top=args$top_n, dir=args$dir, args$k_reg_genes, args$n_regulons)
    regulons <- lapply(vregulons, function(vr) names(vr$tfmode))
    
    log_message(sprintf("GENIE3 regulons created in %s", format_runtime(t0)), logfile, 2)
    t0 <- Sys.time()
    
    modules_exist <- TRUE
    if (length(names(regulons)) == 0) {
      log_message("With the chosen parameters no regulons could be created.", logfile, 2)
      modules_exist <- FALSE
    } else {
      write_json(regulons, path = file.path(args$dir, paste0("genie_genesets.json")), pretty = TRUE, auto_unbox = TRUE)
      
      log_message(sprintf("GENIE3 regulons written in %s", format_runtime(t0)), logfile, 2)
    }
    
    if (modules_exist & (args$oe | args$gsva | args$ssgsea | args$aucell)) {
      
      if (args$oe) {
        log_message("Calculating SPONGeffects-OE scores (GENIE3) ...", logfile)
        t0 <- Sys.time()
        
        oe_scores <- calc_spongeffects_oe(expr, regulons, n_cores=args$n_cores)
        
        log_message(sprintf("SPONGeffects-OE scores computed in %s", format_runtime(t0)), logfile, 2)
        t0 <- Sys.time()
        
        write.csv(oe_scores, file = file.path(args$dir, paste0("spongeffects_OE_scores_genie3.csv")), row.names = TRUE)
        
        log_message(sprintf("SPONGeffects-OE scores written in %s", format_runtime(t0)), logfile, 2)
      }
      
      if (args$gsva) {
        log_message("Calculating SPONGeffects-GSVA scores (GENIE3) ...", logfile)
        t0 <- Sys.time()
        
        gsva_scores <- calc_spongeffects_gsva(expr, regulons, n_cores=args$n_cores)
        
        log_message(sprintf("SPONGeffects-GSVA scores computed in %s", format_runtime(t0)), logfile, 2)
        t0 <- Sys.time()
        
        write.csv(gsva_scores, file = file.path(args$dir, paste0("spongeffects_GSVA_scores_genie3.csv")), row.names = TRUE)
        
        log_message(sprintf("SPONGeffects-GSVA scores written in %s", format_runtime(t0)), logfile, 2)
      }
      
      if (args$ssgsea) {
        log_message("Calculating SPONGeffects-ssGSEA scores (GENIE3) ...", logfile)
        t0 <- Sys.time()
        
        ssgsea_scores <- calc_spongeffects_ssgsea(expr, regulons, n_cores=args$n_cores)
        
        log_message(sprintf("SPONGeffects-ssGSEA scores computed in %s", format_runtime(t0)), logfile, 2)
        t0 <- Sys.time()
        
        write.csv(ssgsea_scores, file = file.path(args$dir, paste0("spongeffects_ssGSEA_scores_genie3.csv")), row.names = TRUE)
        
        log_message(sprintf("SPONGeffects-ssGSEA scores written in %s", format_runtime(t0)), logfile, 2)
      }
      
      # calculate AUCell scores
      if (args$aucell) {
        log_message("Calculating AUCell scores (GENIE3) ...", logfile)
        t0 <- Sys.time()
        
        aucell_scores <- calc_aucell_score(expr, regulons, args$n_cores)
        
        log_message(sprintf("AUCell scores computed in %s", format_runtime(t0)), logfile, 2)
        t0 <- Sys.time()
        
        write.csv(aucell_scores, file = file.path(args$dir, paste0("aucell_scores_genie3.csv")), row.names = TRUE)
        
        log_message(sprintf("AUCell scores written in %s", format_runtime(t0)), logfile, 2)
      }
    }
    
    # calculate VIPER scores
    if (args$viper & modules_exist) {
      log_message("Calculating VIPER scores ...", logfile)
      t0 <- Sys.time()
      
      vpres <- as.data.frame(viper(expr, vregulons, verbose = FALSE))
      
      log_message(sprintf("VIPER scores computed in %s", format_runtime(t0)), logfile, 2)
      t0 <- Sys.time()
      
      write.csv(vpres, file = file.path(args$dir, paste0("viper_scores_genie3.csv")), row.names = TRUE)
      
      log_message(sprintf("VIPER scores written in %s", format_runtime(t0)), logfile, 2)
    }
  }
  
  log_message(paste0("Pipeline finished at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), logfile)
  
}


# only run main() if the script is executed directly, not sourced
if (identical(environment(), globalenv()) && !length(sys.frames())) {
  main()
}
