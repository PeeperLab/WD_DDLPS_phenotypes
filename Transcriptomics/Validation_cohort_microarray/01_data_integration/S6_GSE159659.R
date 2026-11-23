### This script is used to deal with Clariom_S_Human chip data–GSE159659.
## Extract normalized and log2 transformed data from metadata.
library(tibble)
library(dplyr)
library(GEOquery)

## ---- Working directory ----------------------------------------------------
## Before running, set the working directory to the folder holding the
## raw microarray input data (expected folder layout: see README), e.g.
## setwd("/path/to/microarray_raw_data")

out_dir <- "data/rds"
if(!dir.exists(out_dir)){
  dir.create(out_dir, recursive = TRUE) }
  
  
# change the GSE ID if necessary
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 10)
gse <- getGEO("GSE159659", GSEMatrix = TRUE)

# get the ExpressionSet
eset <- gse[[1]]
# extract feature data
feature_data <- fData(eset)

# extract the gene name in "SPOT_ID...10" 
extracted_part <- sub(".*//.*?//(.*?)//.*", "\\1", feature_data$SPOT_ID.1)
feature_data$`Gene Symbol` <- sub(".*\\(([^()]*)\\)[^()]*$", "\\1", extracted_part)
feature_data$`Gene Symbol`[!grepl("\\(.*\\)", extracted_part)] <- "NA"

filtered_feature_data <- feature_data[,c("ID","Gene Symbol")]

# get the expression matrix from the metadata
gene_expression_data <- as.data.frame(exprs(eset))

gene_expression_data <- gene_expression_data %>%
  rownames_to_column(var = "probe")

merged_data <- left_join(filtered_feature_data,gene_expression_data, by = c("ID" = "probe"))
colnames(merged_data)[2] <- "SYMBOL"

saveRDS(merged_data, file.path(out_dir,"gene_exprs_GSE159659.rds"))

