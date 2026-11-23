# This script performs batch correction of multiple proteomic batches using ComBat

library(openxlsx)
library(dplyr)
library(sva)

set.seed(1507)
## ---- Working directory ----------------------------------------------------
## Before running, set the working directory to the "Proteomics" folder of
## this repository (the folder that holds Inputs/ and Outputs/), e.g.
## setwd("/path/to/WD_DDLPS_phenotypes/Proteomics")

protein_data <- read.xlsx("Inputs/Proteomics_Integrated.xlsx")
metadata <- read.xlsx("Inputs/Proteomics_metadata.xlsx")

output_dir <- "Inputs"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Put ID as rownames
rownames(protein_data) <- protein_data[[1]]
protein_data[[1]] <- NULL
rownames(metadata) <- metadata$ID

# Match samples by ID (not by position). Only samples with metadata (batch)
# can be batch-corrected. Sample 48_F (healthy fat, batch 1172) is present in
# the proteomic data but has no entry in Proteomic_metadata.xlsx and is
# therefore excluded (189 -> 188 samples). All downstream proteomic analyses
# and figures use these 188 samples.
common_ids <- sort(intersect(rownames(protein_data), metadata$ID))
dropped <- setdiff(rownames(protein_data), common_ids)
if (length(dropped) > 0) {
  message("Samples without metadata, excluded from batch correction: ",
          paste(dropped, collapse = ", "))
}
protein_data <- protein_data[common_ids, ]
metadata <- metadata[common_ids, ]
stopifnot(identical(rownames(protein_data), metadata$ID))

######################## Batch correction ################################
# transform the data
scaled_data <- t(log2(protein_data + 1))   # proteins x samples
batch_info <- as.character(metadata$batch)

# Perform ComBat batch correction
corrected_data <- ComBat(dat = scaled_data, batch = batch_info, mod = NULL)
corrected_data <- t(corrected_data)

write.csv(corrected_data, file = file.path(output_dir, "Proteomic_data_batch_corrected.csv"))

