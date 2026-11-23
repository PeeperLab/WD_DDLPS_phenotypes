# This script performs differential expression analysis using Limma between all clusters identified in the consensus clustering.

library(readxl)
library(limma)
library(openxlsx)  # read.xlsx()

set.seed(1507)
## ---- Working directory ----------------------------------------------------
## Before running, set the working directory to the "Proteomics" folder of
## this repository (the folder that holds Inputs/ and Outputs/), e.g.
## setwd("/path/to/WD_DDLPS_phenotypes/Proteomics")

proteomic_data <- read.csv("Inputs/Proteomic_data_batch_corrected.csv", row.names = 1)
proteomic_metadata <- read.xlsx("Inputs/Proteomics_metadata.xlsx")

output_dir <- "Outputs/Files/Limma"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Put ID as rownames and order
rownames(proteomic_metadata) <- proteomic_metadata[[1]]   # use first column as rownames

proteomic_metadata <- proteomic_metadata[
  proteomic_metadata$class == "DD" | proteomic_metadata$class == "WD", ]

proteomic_data <- proteomic_data[rownames(proteomic_data) %in% rownames(proteomic_metadata), ]
proteomic_data <- t(proteomic_data)

clusters <- read.csv("Outputs/Files/Consensus_clustering/Consensus_clustering.csv", row.names = 1)

# Create a dataframe from the DD consensus clusters
cluster_df <- data.frame(
  Cluster = clusters$consensus_clusters,  # Values from your vector
  row.names = rownames(proteomic_metadata)  # Set rownames
)

# Initialize lists to store results
ranked_gene_lists <- list()
de_proteins_list <- list()

# Loop through each cluster
for (cluster in sort(unique(cluster_df$Cluster))) {
  
  # Define binary contrast: 1 = in cluster, 0 = all others
  cluster_vs_all <- ifelse(cluster_df$Cluster == cluster, 1, 0)
  
  # Create design matrix
  design <- model.matrix(~ cluster_vs_all)
  
  # Fit linear model
  fit <- lmFit(proteomic_data, design)
  fit <- eBayes(fit)
  
  # Extract ranked gene list (t-statistics)
  ranked_gene_lists[[paste0("Cluster_", cluster)]] <-
    setNames(fit$t[, 2], rownames(fit))
  
  # Extract ALL DE results
  de_results <- topTable(
    fit,
    coef = 2,
    adjust.method = "fdr",
    number = Inf
  )
  
  de_proteins_list[[paste0("Cluster_", cluster)]] <- de_results
  
  # Define output file
  out_file <- file.path(
    output_dir,
    paste0("Cluster_", cluster, "_DE_proteins_6cl.csv")
  )
  
  # Save CSV
  write.csv(
    de_results,
    file = out_file,
    row.names = TRUE
  )
}
