###This SCRIPT is used to change the microarray data to gene expression level.

library(affy)
library(GEOquery)
library(tidyverse)
library(hgu133plus2.db)
library(hgu133a.db)
library(tibble)

# set working directory
## ---- Working directory ----------------------------------------------------
## Before running, set the working directory to the folder holding the
## raw microarray input data (expected folder layout: see README), e.g.
## setwd("/path/to/microarray_raw_data")

gse_ids <- c("GSE21050", "GSE21122", "GSE30929", "GSE71118")
lps_root <- "data/lps_data"

out_dir <- "data/rds"

for (gse_id in gse_ids) {
  #read .CEL files
  cel_files_path <- file.path(lps_root,paste0(gse_id,"_LPS"))
  cel_files <- list.files(path = cel_files_path, pattern = "*.CEL.gz", full.names = TRUE)
  
  # Read the .CEL files into an AffyBatch object
  raw_data <- ReadAffy(filenames = cel_files)
  
  # remove the .CEL.gz in the sample name
  sample_names <- gsub("\\.CEL\\.gz", "", sampleNames(raw_data))
  sampleNames(raw_data) <- sample_names
  
  # correct background without prior normalization
  normalized_data <- affy::rma(raw_data)
  
  # obtain the gene expression matrix
  expression_data <- as.data.frame(exprs(normalized_data))
  expression_data <- expression_data %>%
    rownames_to_column(var = "probe")
  
  probe_ids <- expression_data$probe
  
  if (gse_id %in% c("GSE21050","GSE71118")) {
    gene_symbols <- mapIds(hgu133plus2.db,  #change this chip information based the chip they use
                           keys = probe_ids,  
                           column = "SYMBOL",  
                           keytype = "PROBEID", 
                           multiVals = "first") 
  } else {
    gene_symbols <- mapIds(hgu133a.db,  #change this chip information based the chip they use
                           keys = probe_ids,  
                           column = "SYMBOL",  
                           keytype = "PROBEID", 
                           multiVals = "first") 
    
      }
  gene_symbols <- as.data.frame(gene_symbols)
  
  gene_symbols <- gene_symbols %>%
    rownames_to_column(var="probe")
  
  merged_data <- left_join( gene_symbols,expression_data, by = "probe")
  
  ## change the file name based on the input file
  saveRDS(merged_data, file.path(out_dir,paste0("gene_exprs_", gse_id,".rds")))
 }