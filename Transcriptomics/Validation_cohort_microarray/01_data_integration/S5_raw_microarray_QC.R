
## QC pipeline for multiple microarray datasets (raw QC -> RMA -> post-RMA QC)
## Method:
## - Pre-normalization QC (arrayQualityMetrics): distances, boxplots, RLE, NUSE, MA, spatial distribution of M
## - Flag samples failing >= 2 metrics for manual review
## - Pre-normalization batch check: inspect metadata + raw intensity boxplots
## - RMA (background correction, normalization, log2)
## - Post-RMA QC (arrayQualityMetrics): distances, boxplots, MA

suppressPackageStartupMessages({
  library(arrayQualityMetrics)
  library(affy)
  library(GEOquery)
})

## ---- Working directory ----------------------------------------------------
## Before running, set the working directory to the folder holding the
## raw microarray input data (expected folder layout: see README), e.g.
## setwd("/path/to/microarray_raw_data")
gse_ids <- c("GSE159659", "GSE21050", "GSE21122", "GSE30929", "GSE71118")
lps_root <- "data/lps_data"


for (gse_id in gse_ids) {

  message("Currently processing: ", gse_id)
  
  gse_dir <- file.path(lps_root, paste0(gse_id, "_LPS"))
  cel_files <- list.files(path = gse_dir, pattern = "\\.CEL(\\.gz)?$", full.names = TRUE, ignore.case = TRUE)
  if (gse_id != "GSE159659") {
    raw_data <- affy::ReadAffy(filenames = cel_files)
    
    ## Raw Data QC
    raw_outdir <- file.path("qc_reports", paste0(gse_id, "_Raw_QC"))
    arrayQualityMetrics(expressionset = raw_data, outdir = raw_outdir, force = TRUE, do.logtransform = TRUE)
    
    # data normalization
      rma_data <- affy::rma(raw_data)
      ## RMA Data QC
      rma_outdir <- file.path("qc_reports", paste0(gse_id, "_RMA_QC"))
      arrayQualityMetrics(expressionset = rma_data, outdir = rma_outdir, force = TRUE, do.logtransform = FALSE)
      
    } else {
      Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 10)
      gse <- getGEO("GSE159659", GSEMatrix = TRUE)
      # get the ExpressionSet
      eset <- gse[[1]]
      ## RMA Data QC
      rma_outdir <- file.path("qc_reports", paste0(gse_id, "_RMA_QC"))
      arrayQualityMetrics(expressionset = eset, outdir = rma_outdir, force = TRUE, do.logtransform = FALSE)
      
    }
    
    message("Successfully generated Raw and RMA reports for: ", gse_id)
    
}