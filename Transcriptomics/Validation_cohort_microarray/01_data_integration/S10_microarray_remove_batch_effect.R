## This script is used to remove the batch effect.

library(sva)
library(dplyr)
library(readr)
library(tibble)

# set working directory
## ---- Working directory ----------------------------------------------------
## Before running, set the working directory to the folder holding the
## raw microarray input data (expected folder layout: see README), e.g.
## setwd("/path/to/microarray_raw_data")

file_path1 <- "meta/metadata/Metadata_validation_all.csv" 
file_path2 <-  "data/rds/integrated_qc_exprs_microarray.rds"

out_dir <- "data/rds"
if(!dir.exists(out_dir)){
  dir.create(out_dir, recursive = TRUE) }
## metadata
metadata <- read_csv(file_path1)
metadata$Dataset <- metadata$dataset
# expression data
raw_exprs_data <- readRDS(file_path2)
exprs_data <- as.data.frame(raw_exprs_data)
rownames(exprs_data) <- exprs_data$SYMBOL
exprs_data <- exprs_data[,-1]
exprs_data <- exprs_data[,metadata$sampleID]
#all(metadata$sampleID == colnames(exprs_data))

##### remove batch effect
batch_vector <- as.factor(metadata$dataset)
#  remove batch effect
combat_data <- ComBat(exprs_data, batch = batch_vector, par.prior=TRUE, prior.plots=FALSE)
# data.frame
combat_data_df <- as.data.frame(combat_data)
saveRDS(combat_data_df, file = file.path(out_dir,"Batch_corrected_data_validation_all.rds"))

