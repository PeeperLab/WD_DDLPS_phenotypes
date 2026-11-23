### This script combines all expression data and harmonizes gene symbols.

library(dplyr)
library(readr)

# set working directory
## ---- Working directory ----------------------------------------------------
## Before running, set the working directory to the folder holding the
## raw microarray input data (expected folder layout: see README), e.g.
## setwd("/path/to/microarray_raw_data")

input_dir <- "data/rds"

# select raw expression rds files
rds_files <- list.files(
  path = input_dir,
  pattern = "^gene_exprs_.*\\.rds$",
  full.names = TRUE
)


# low-quality samples to remove at the end
low_qc_sample <- c("GSM528276","GSM766533","GSM1827111","GSM4837054","GSM4837069")

processed_list <- list()

for (f in rds_files) {
  suffix <- sub("^gene_exprs_|\\.rds$", "", basename(f))
  
  expression_data <- readRDS(f)
  
  # make sure it's a data.frame
  expression_data <- as.data.frame(expression_data)
  expression_data <- expression_data[ ,-1]
  colnames(expression_data)[1] <- "SYMBOL"
  # drop NA SYMBOL
  expression_data <- expression_data %>%
    filter(!is.na(SYMBOL))
  
  # if num !=0,delete the rows with GeneSymbol is NA
  expression_data <- expression_data[!is.na(expression_data[,1]), ]
  
 # expression_data <- expression_data[,-1]
  
  expression_data_unique <- expression_data %>%
    group_by(SYMBOL) %>%
    summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))
  
  ## save per-dataset file
  #out_file <- file.path(input_dir, paste0("unique_gene_expr_", suffix))
  # saveRDS(expression_data_unique, out_file)
  # message("  Saved unique file: ", out_file, "  (dim=", paste(dim(expression_data_unique), collapse="x"), ")")

  processed_list[[suffix]] <- expression_data_unique
}

# merge all processed data by SYMBOL (intersection of genes across datasets)
merged_exprs <- Reduce(function(x, y) inner_join(x, y, by = "SYMBOL"), processed_list)

# remove low-quality samples (columns), keep SYMBOL
merged_exprs_qc <- merged_exprs[ ,!colnames(merged_exprs) %in% low_qc_sample]

# save integrated result in input_dir
final_file <- file.path(input_dir, "integrated_qc_exprs_microarray.rds")
saveRDS(merged_exprs_qc, final_file)


