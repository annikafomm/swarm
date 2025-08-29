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

source("./calc_R_scores/utils.R")

#  Rscript calc_scores.R   --dir ../datasets_scores/GSM6592049_M2_prepro   --dir ../datasets_scores/GSM6592049_M2_scores   --sponge_network ../networks/SPONGE/breast_invasive_carcinoma/breast_invasive_carcinoma_interactionNetwork.csv   --sponge_analysis ../networks/SPONGE/breast_invasive_carcinoma/breast_invasive_carcinoma_networkAnalysis.csv --genie_network ../networks/GENIE3/BRCA/genie3_BRCA_mrn.top_100k.csv --ensembl_col ensemble_id  --aucell --gsva --ssgsea --viper

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

compute_network_scores <- function(description, args, logfile) {
  # Load the data
  if (description == "tg") {
    log_message("Preparing score calculation for the Tangram output ...", logfile)
  } else if (description == "st") {
    log_message("Preparing score calculation for the Spatial data ...", logfile)
  }
  
  # expression matrix
  log_message("Loading expression matrix (and additional files) ...", logfile, 2)
  t0 <- Sys.time()
  if (description == "tg") {
    dir_expr <- file.path(args$dir, "expr_info_tg")
  } else if (description == "st"){
    dir_expr <- file.path(args$dir, "expr_info_st")
  }

  expr <- readMM(file.path(dir_expr, "expr.mtx"))
  rownames(expr) <- fread(file.path(dir_expr, "cells.txt"), header = FALSE)$V1
  var_df <- fread(file.path(dir_expr, "var.csv"))
  expr <- t(as.matrix(expr))

  if (nrow(expr) != nrow(var_df)) {
    log_message("The expression matrix and the var file don't belong together.", logfile, 2)
  } else {
    log_message(sprintf("Expression matrix loaded in %s", format_runtime(t0)), logfile, 4)

    # SPONGE networks
    if (!is.null(args$sponge_network) & !is.null(args$sponge_analysis) & (args$aucell | args$gsva | args$ssgsea)) {
      log_message("Loading SPONGE network and create SPONGE modules ...", logfile, 2)
      t0 <- Sys.time()
      # load files and calculate SPONGE modules
      ceRNA_interactions <- fread(args$sponge_network)
      ceRNA_centralities <- fread(args$sponge_analysis)
      log_message(sprintf("SPONGE network loaded in %s", format_runtime(t0)), logfile, 4)

      interaction_cols <- c("geneA", "geneB", "df", "cor", "pcor", "mscor", "p.val", "p.adj")
      centrality_cols <- c("gene", "degree", "eigenvector", "betweenness", "page_rank")
      right_order <- all(interaction_cols %in% colnames(ceRNA_interactions) & all(centrality_cols %in% colnames(ceRNA_centralities)))
      wrong_order <- all(centrality_cols %in% colnames(ceRNA_interactions) & all(interaction_cols %in% colnames(ceRNA_centralities)))
      if (!(right_order | wrong_order)) {
        log_message("At least one of the sponge network files has an unsupported file format.", logfile)
      } else if (! args$ensembl_col %in% colnames(var_df)){
        log_message(paste0("'", args$ensembl_col, "' is not a column in adata.vars. Please provide a valid ensembl column."), logfile)
      } else if (! args$feature_col %in% colnames(var_df)){
        log_message(paste0("'", args$feature_col, "' is not a column in adata.vars. Please provide a valid feature column."), logfile)
        # TODO: check RNAs
      } else {
        if (wrong_order) {
          centralities <- ceRNA_centralities
          ceRNA_centralities <- ceRNA_interactions
          ceRNA_interactions <- centralities
        }

        rownames(expr) <- var_df[,get(args$ensembl_col)]

        t0 <- Sys.time()
        sponge_modules <- create_Sponge_modules(var_df, ceRNA_interactions, ceRNA_centralities, mscor=args$mscor, padj=args$padj,
                                          feature_col=args$feature_col, RNAs=args$RNA_types, ensembl_col=args$ensembl_col, max_modules=args$max_modules, n_cores=args$n_cores)
        log_message(sprintf("SPONGE modules created in %s", format_runtime(t0)), logfile, 4)


        if (length(names(sponge_modules)) == 0) {
          log_message("With the chosen parameters no SPONGE modules could be created.", logfile)
        } else {

          if (description == "tg") {
            outdir_scores <- file.path(args$dir, "Rscores_tg")
          } else if (description == "st") {
            outdir_scores <- file.path(args$dir, "Rscores_st")
          }
          if (!dir.exists(outdir_scores)) {
            dir.create(outdir_scores, recursive = TRUE, showWarnings = FALSE)
          }

          t0 <- Sys.time()
          write_json(sponge_modules, path = file.path(outdir_scores, "sponge_genesets.json"), pretty = TRUE, auto_unbox = TRUE)
          log_message(sprintf("SPONGE modules written in %s", format_runtime(t0)), logfile, 4)

          t0 <- Sys.time()
          filt_sponge_network <- filter_network(ceRNA_interactions, sponge_modules, "Sponge")
          fwrite(filt_sponge_network, file.path(args$dir, paste0("sponge_network_filtered_", description, ".csv")))
          log_message(sprintf("SPONGE network filtered and written in %s", format_runtime(t0)), logfile, 4)

          log_message("Calculating scores with the SPONGE network ...", logfile, 2)

          # calculate SPONGeffects scores

          if (args$gsva) {
            t0 <- Sys.time()
            gsva_scores <- calc_spongeffects_gsva(expr, sponge_modules, n_cores=args$n_cores)
            log_message(sprintf("SPONGeffects-GSVA scores computed in %s", format_runtime(t0)), logfile, 4)

            t0 <- Sys.time()
            write.csv(gsva_scores, file = file.path(outdir_scores, "spongeffects_GSVA_scores_sponge.csv"), row.names = TRUE)
            log_message(sprintf("SPONGeffects-GSVA scores written in %s", format_runtime(t0)), logfile, 4)
          }

          if (args$ssgsea) {
            t0 <- Sys.time()
            ssgsea_scores <- calc_spongeffects_ssgsea(expr, sponge_modules, n_cores=args$n_cores)
            log_message(sprintf("SPONGeffects-ssGSEA scores computed in %s", format_runtime(t0)), logfile, 4)
            
            t0 <- Sys.time()
            write.csv(ssgsea_scores, file = file.path(outdir_scores, "spongeffects_ssGSEA_scores_sponge.csv"), row.names = TRUE)
            log_message(sprintf("SPONGeffects-ssGSEA scores written in %s", format_runtime(t0)), logfile, 4)
          }

          # calculate AUCell scores
          if (args$aucell) {
            t0 <- Sys.time()
            aucell_scores <- calc_aucell_score(expr, sponge_modules, args$n_cores)
            log_message(sprintf("AUCell scores computed in %s", format_runtime(t0)), logfile, 4)
            
            t0 <- Sys.time()
            write.csv(aucell_scores, file = file.path(outdir_scores, "aucell_scores_sponge.csv"), row.names = TRUE)
            log_message(sprintf("AUCell scores written in %s", format_runtime(t0)), logfile, 4)
          } 
        }
      }
    }

    # GENIE3 networks

    if (! is.null(args$genie_network)) {
      log_message("Loading GENIE3 network and create regulons ...", logfile, 2)

      rownames(expr) <- var_df[,V1]

      t0 <- Sys.time()
      # load network file and calculate regulons
      regulon_df <- fread(args$genie_network)
      log_message(sprintf("GENIE3 network loaded in %s", format_runtime(t0)), logfile, 4)

      if (description == "tg") {
        outdir_scores <- file.path(args$dir, "Rscores_tg")
      } else if (description == "st") {
        outdir_scores <- file.path(args$dir, "Rscores_st")
      }
      if (!dir.exists(outdir_scores)) {
        dir.create(outdir_scores, recursive = TRUE, showWarnings = FALSE)
      }
      
      t0 <- Sys.time()
      vregulons <- create_Genie_modules(regulon_df, expr=expr, top=args$top_n, dir=outdir_scores, args$k_reg_genes, args$n_regulons)
      regulons <- lapply(vregulons, function(vr) names(vr$tfmode))
      log_message(sprintf("GENIE3 regulons created in %s", format_runtime(t0)), logfile, 4)

      if (length(names(regulons)) == 0) {
        log_message("With the chosen parameters no regulons could be created.", logfile)
      } else {
        t0 <- Sys.time()
        write_json(regulons, path = file.path(outdir_scores, "genie_genesets.json"), pretty = TRUE, auto_unbox = TRUE)
        log_message(sprintf("GENIE3 regulons written in %s", format_runtime(t0)), logfile, 4)
        
        t0 <- Sys.time()
        filt_genie_network <- filter_network(regulon_df, regulons, "Genie")
        fwrite(filt_genie_network, file.path(args$dir, paste0("genie_network_filtered_", description, ".csv")))
        log_message(sprintf("GENIE3 network filtered and written in %s", format_runtime(t0)), logfile, 4)

        log_message("Calculating scores with the GENIE network ...", logfile, 2)

        # calculate SPONGeffects scores

        if (args$gsva) {
          t0 <- Sys.time()
          gsva_scores <- calc_spongeffects_gsva(expr, regulons, n_cores=args$n_cores)
          log_message(sprintf("SPONGeffects-GSVA scores computed in %s", format_runtime(t0)), logfile, 4)
          
          t0 <- Sys.time()
          write.csv(gsva_scores, file = file.path(outdir_scores, "spongeffects_GSVA_scores_genie3.csv"), row.names = TRUE)
          log_message(sprintf("SPONGeffects-GSVA scores written in %s", format_runtime(t0)), logfile, 4)
        }

        if (args$ssgsea) {
          t0 <- Sys.time()
          ssgsea_scores <- calc_spongeffects_ssgsea(expr, regulons, n_cores=args$n_cores)
          log_message(sprintf("SPONGeffects-ssGSEA scores computed in %s", format_runtime(t0)), logfile, 4)
          
          t0 <- Sys.time()
          write.csv(ssgsea_scores, file = file.path(outdir_scores, "spongeffects_ssGSEA_scores_genie3.csv"), row.names = TRUE)
          log_message(sprintf("SPONGeffects-ssGSEA scores written in %s", format_runtime(t0)), logfile, 4)
        }

        # calculate AUCell scores
        if (args$aucell) {
          t0 <- Sys.time()
          aucell_scores <- calc_aucell_score(expr, regulons, args$n_cores)
          log_message(sprintf("AUCell scores computed in %s", format_runtime(t0)), logfile, 4)
          
          t0 <- Sys.time()
          write.csv(aucell_scores, file = file.path(outdir_scores, "aucell_scores_genie3.csv"), row.names = TRUE)
          log_message(sprintf("AUCell scores written in %s", format_runtime(t0)), logfile, 4)
        }

        # calculate VIPER scores
        if (args$viper) {
          t0 <- Sys.time()
          vpres <- as.data.frame(viper(expr, vregulons, verbose = FALSE))
          log_message(sprintf("VIPER scores computed in %s", format_runtime(t0)), logfile, 4)
          
          t0 <- Sys.time()
          write.csv(vpres, file = file.path(outdir_scores, "viper_scores_genie3.csv"), row.names = TRUE)
          log_message(sprintf("VIPER scores written in %s", format_runtime(t0)), logfile, 4)
        }
      }
    }
  }
}



main <- function() {
  option_list <- list(
    # dir with files
    make_option("--dir", type="character", help="Dir containing all input files (expr.mtx, var.tsv, (cells.txt))"),
    make_option("--log", type="character", help="Path to the log file", metavar="file"),

    make_option("--tangram", action="store_true", default=FALSE, help="Was tangram used?"),

    # network files
    make_option("--sponge_network", type="character", help="Path to Sponge interactionNetwork file", metavar="file"),
    make_option("--sponge_analysis", type="character", help="Path to Sponge networkAnalysis file", metavar="file"),
    make_option("--genie_network", type="character", help="Path to GENIE3 network file", metavar="file"),

    # score flags
    make_option("--aucell", action="store_true", default=FALSE, help="Calculate AUCell score"),
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

  # required flags
  if (is.null(args$log)) {
    stop("Error: --log is required.")
  }
  if (is.null(args$dir)) {
    stop("Error: --dir is required.")
  }
  if (! all(c("expr.mtx", "var.csv", "cells.txt") %in% list.files(file.path(args$dir, "expr_info_st")))) {
    stop(paste0("Error: The directory ", args$dir, " does not contain the dir 'expr_info_st' with the files (expr.mtx, var.csv, cells.txt)"))
  }
  if (args$tangram & ! all(c("expr.mtx", "var.csv", "cells.txt") %in% list.files(file.path(args$dir, "expr_info_tg")))) {
    stop(paste0("Error: The directory ", args$dir, " does not contain the dir 'expr_info_tg' with the files (expr.mtx, var.csv, cells.txt)"))
  }
  if ((is.null(args$sponge_network) | is.null(args$sponge_analysis)) & is.null(args$genie_network)) {
    stop("Error: A SPONGE (--sponge_network and --sponge_analysis) or a GENIE3 (--genie_network) network is required.")
  }

  # Prepare log file
  logfile <- args$log
  log_message(paste0("R score pipeline started at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), logfile)

  if (args$tangram) {
    compute_network_scores("tg", args, logfile)
  }
  compute_network_scores("st", args, logfile)
  
  log_message(paste0("R score pipeline finished at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n"), logfile)
}


# only run main() if the script is executed directly, not sourced
if (identical(environment(), globalenv()) && !length(sys.frames())) {
  main()
}
