### This script selects WDLPS/DDLPS samples from the integrated (RMA-normalised) microarray data
### and corrects the batch (dataset) effect with ComBat.
set.seed(123)
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
out_dir1 <- "meta/metadata"
out_dir2 <- "data/rds"
for (d in c(out_dir1, out_dir2)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

### metadata
metadata <- read_csv(file_path1)

metadata <- metadata[metadata$subtype %in% c("WDLPS","DDLPS"),]
write_csv(metadata,file = file.path(out_dir1,"Metadata_validation_WD_DD.csv"))
### couunt data
count_data <- readRDS(file_path2)
count_data_unique <- as.data.frame(count_data)
# change the gene_symbol column to row names
rownames(count_data_unique) <- count_data_unique$SYMBOL
count_data_unique$SYMBOL <- NULL
count_data_unique <- count_data_unique[,metadata$sampleID]
count_data <- count_data_unique

#####remove batch effect
batch_vector <- as.factor(metadata$dataset)
#  remove batch effect
combat_data <- ComBat(count_data, batch = batch_vector, par.prior=TRUE, prior.plots=FALSE)
# data.frame
combat_data_df <- as.data.frame(combat_data)
saveRDS(combat_data_df, file = file.path(out_dir2,"Batch_corrected_data_validation_WD_DD.rds"))
