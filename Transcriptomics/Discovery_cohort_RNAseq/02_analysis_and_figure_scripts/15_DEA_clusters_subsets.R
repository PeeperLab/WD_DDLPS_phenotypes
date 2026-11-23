# This script performs differential expression analysis via Limma between specific clusters

library(readr)
library(dplyr)
library(limma)

# set working directory
set.seed(1507)
## ---- Working directory ----------------------------------------------------
## Before running, set the working directory to the "Transcriptomics/Discovery_cohort_RNAseq" folder of
## this repository (the folder that holds Inputs/ and Outputs/), e.g.
## setwd("/path/to/WD_DDLPS_phenotypes/Transcriptomics/Discovery_cohort_RNAseq")

output_dir <- "Outputs/Files/Limma/Cluster_1_vs_2"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Load data
metadata <- read_csv("Inputs/Metadata_discovery_WD_DD.csv")
exprs_data    <- readRDS("Inputs/Batch_corrected_data_discovery_WD_DD.rds")
clusters6 <- read.csv("Outputs/Files/Consensus_clustering/Consensus_clustering.csv", row.names = 1)

metadata <- metadata %>%
  left_join(clusters6, by = "sampleID")
colnames(metadata)[ncol(metadata)] <- "clusters6"

# Select clusters of interest
clusters_of_interest <- c("1", "2")
metadata_sub <- metadata %>%
  filter(clusters6 %in% clusters_of_interest)

exprs_sub <- exprs_data[, metadata_sub$sampleID]

# DEA: cluster vs rest (within subset)
dea_results_list <- list()

for (cluster in sort(unique(metadata_sub$clusters6))) {
  
  cat("Performing DEA for cluster", cluster, "...\n")
  
  cluster_vs_rest <- ifelse(metadata_sub$clusters6 == cluster, 1, 0)
  design <- model.matrix(~ cluster_vs_rest)
  
  fit <- lmFit(exprs_sub, design)
  fit <- eBayes(fit)
  
  de_results <- topTable(
    fit,
    coef = 2,
    adjust.method = "fdr",
    number = Inf
  )
  
  # Save DEA results CSV
  out_csv <- file.path(
    output_dir,
    paste0("Cluster_", cluster, "_vs_", paste(setdiff(clusters_of_interest, cluster), collapse = "_"), "_DEGs.csv")
  )
  write.csv(de_results, out_csv, row.names = TRUE)
  
  dea_results_list[[paste0("Cluster", cluster)]] <- de_results
}
