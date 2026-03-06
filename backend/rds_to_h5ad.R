#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(SeuratObject)
  library(SeuratDisk)
  library(optparse)
})


main <- function() {

  option_list <- list(
    make_option(
      c("--rds_path"),
      type = "character",
      help = "Path to Seurat .rds file (required)."
    ),
    make_option(
      c("--assay"),
      type = "character",
      default = "RNA",
      help = "Assay name to export (default: %default)."
    ),
    make_option(
      c("--h5ad_path"),
      type = "character",
      default = "",
      help = "Path to which h5ad file will be written (default: same as --rds_path with .h5ad suffix)."
    )
  )

  opt_parser <- OptionParser(option_list = option_list)
  opt <- parse_args(opt_parser)

  cat("Package versions:\n")
  cat("  SeuratObject:", as.character(packageVersion("SeuratObject")), "\n")
  cat("  Seurat:", as.character(packageVersion("Seurat")), "\n")
  cat("  SeuratDisk:", as.character(packageVersion("SeuratDisk")), "\n")
  cat("  Signac:", as.character(packageVersion("Signac")), "\n")

  cat("SeuratObject built:", packageDescription("SeuratObject")$Built, "\n")
  cat("SeuratObject:", as.character(packageVersion("SeuratObject")), "\n")
  cat("SeuratObject path:", find.package("SeuratObject"), "\n")
  cat("Seurat:", as.character(packageVersion("Seurat")), "\n")
  cat("Seurat path:", find.package("Seurat"), "\n")

  cat("SeuratDisk path:", find.package("SeuratDisk"), "\n")
  cat("Signac path:", find.package("Signac"), "\n")


  # ---- Basic argument checks ----
  if (is.null(opt$rds_path)) {
    stop("Missing required argument: --rds_path\n", call. = FALSE)
  }

  rds_path   <- opt$rds_path
  assay_name <- opt$assay
  h5ad_path  <- opt$h5ad_path

  cat("Reading Seurat object from:", rds_path, "\n")
  data_seurat <- readRDS(rds_path)

  if (!assay_name %in% names(data_seurat@assays)) {
    stop(
      sprintf(
        "Assay '%s' not found in Seurat object. Available assays: %s",
        assay_name,
        paste(names(data_seurat@assays), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  # ---- Build output file names ----
  h5seurat_file <- sub("\\.rds$", ".h5Seurat", rds_path)
  if (h5ad_path==""){
    h5ad_file <- sub("\\.h5Seurat$", paste0("_", assay_name, ".h5ad"), h5seurat_file)
  }
  else{
    h5ad_file <- h5ad_path
  }

  cat("h5Seurat file will be:", h5seurat_file, "\n")
  cat("h5ad file will be:    ", h5ad_file, "\n")

  # # ---- Ensure assay is standard Assay class ----
  # cat("Coercing assay", assay_name, "to 'Assay' class if needed...\n")
  # data_seurat[[assay_name]] <- as(object = data_seurat[[assay_name]], Class = "Assay")
  # ---- Convert assay to V3/4 Assay structure using scCustomize ----
  cat("Converting assay", assay_name, "to V3/4 Assay structure via scCustomize::Convert_Assay()...\n")

  if (!requireNamespace("scCustomize", quietly = TRUE)) {
    stop("Package 'scCustomize' is required for option 2. Install with: install.packages('scCustomize')", call. = FALSE)
  }

  # Convert only the assay you plan to export
  data_seurat <- scCustomize::Convert_Assay(
    seurat_object = data_seurat,
    assay         = assay_name,
    convert_to    = "V3"   # accepted: "V3"/"Assay"/"v3"/"assay"
  )

  # (Optional but often helpful) set the default assay explicitly
  Seurat::DefaultAssay(data_seurat) <- assay_name


  # ---- Save to .h5Seurat ----
  cat("Saving Seurat object as .h5Seurat...\n")
  SeuratDisk::SaveH5Seurat(
    object    = data_seurat,
    filename  = h5seurat_file,
    overwrite = TRUE
  )

  # ---- Convert to .h5ad (for Scanpy) ----
  cat("Converting .h5Seurat to .h5ad (assay:", assay_name, ")...\n")
  SeuratDisk::Convert(
    h5seurat_file,
    dest      = "h5ad",
    assay     = assay_name,
    filename  = h5ad_file,
    overwrite = TRUE
  )

  cat("Done. Wrote h5ad to:", h5ad_file, "\n")
}

if (sys.nframe() == 0) {
  main()
}
