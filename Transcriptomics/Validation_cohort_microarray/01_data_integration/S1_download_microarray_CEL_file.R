###This script is used to download and untar the microarray data.
# You can also download the data from GEO manually !!!
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("GEOquery")

library(GEOquery)

# set working directory
## ---- Working directory ----------------------------------------------------
## Before running, set the working directory to the folder holding the
## raw microarray input data (expected folder layout: see README), e.g.
## setwd("/path/to/microarray_raw_data")


# set the path for the download data
download_path <- "raw_data"
dir.create(download_path, showWarnings = FALSE)

# define the GSE id, change the id for each dataset
gse_ids <- c("GSE21050", "GSE21122", "GSE30929", "GSE71118", "GSE159659")

for (gse_id in gse_ids) {
  # download files
  getGEOSuppFiles(gse_id, baseDir = download_path)
  
  # untar files
  raw_dir <- file.path(download_path, gse_id, paste0(gse_id, "_RAW"))  # read by S4
  if (!dir.exists(raw_dir)) {
    dir.create(raw_dir)
  }
  
  
  tar_file <- file.path(download_path, gse_id, paste0(gse_id, "_RAW.tar"))
  untar(tar_file, exdir = raw_dir)
  
  cat(paste("Complete：", gse_id, "\n"))
}
