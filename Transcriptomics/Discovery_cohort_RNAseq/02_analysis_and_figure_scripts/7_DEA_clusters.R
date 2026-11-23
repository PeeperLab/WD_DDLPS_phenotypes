# This script performs differential gene expression analysis (DEA) for each of the 6 clusters obtained 
# from consensus clustering on the WD/DD subset of the transcriptomic discovery cohort. 

library(readr)
library(dplyr)
library(limma)
library(tibble)

# set working directory
set.seed(1507)
## ---- Working directory ----------------------------------------------------
## Before running, set the working directory to the "Transcriptomics/Discovery_cohort_RNAseq" folder of
## this repository (the folder that holds Inputs/ and Outputs/), e.g.
## setwd("/path/to/WD_DDLPS_phenotypes/Transcriptomics/Discovery_cohort_RNAseq")

output_dir <- "Outputs/Files/Limma"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# load data
metadata <- read_csv("Inputs/Metadata_discovery_WD_DD.csv")
exprs_data <- readRDS("Inputs/Batch_corrected_data_discovery_WD_DD.rds")
k_6clusters <- read.csv("Outputs/Files/Consensus_clustering/Consensus_clustering.csv", row.names = 1)

## add the consensus clustering result to metadata
metadata <- left_join(metadata,k_6clusters, by = "sampleID")
colnames(metadata)[5] <- "clusters6"
metadata <- metadata[match(colnames(exprs_data), metadata$sampleID), ]  # align with expression columns
stopifnot(identical(metadata$sampleID, colnames(exprs_data)))

# Initialize lists to store results
ranked_gene_lists <- list()
de_proteins_list <- list()

# Loop through each cluster
for (cluster in sort(unique(metadata$clusters6))) {
  
  # Define binary contrast: 1 = in cluster, 0 = all others
  cluster_vs_all <- ifelse(metadata$clusters6 == cluster, 1, 0)
  
  # Create design matrix
  design <- model.matrix(~ cluster_vs_all)
  
  # Fit linear model
  fit <- lmFit(exprs_data, design)
  fit <- eBayes(fit)
  
  # Extract ranked gene list (t-statistics)
  ranked_gene_lists[[paste0("Cluster_", cluster)]] <-
    setNames(fit$t[, 2], rownames(fit))
  
  # Extract ALL DE results (full table, no filtering)
  de_results <- topTable(
    fit,
    coef = 2,
    adjust.method = "fdr",
    number = Inf
  )
  
  de_proteins_list[[paste0("Cluster_", cluster)]] <- de_results
  
#  Save the full table to CSV
  write.csv(
    de_results,
    file = file.path(output_dir, paste0("Cluster_", cluster, "_DEGs_all.csv")),
    row.names = TRUE
  )
}

############################## Extracting all DEGs ########################################
# for 6 clusters:
DEA_6cl_1 <- read.csv("Outputs/Files/Limma/Cluster_1_DEGs_all.csv", row.names = 1)
DEA_6cl_2 <- read.csv("Outputs/Files/Limma/Cluster_2_DEGs_all.csv", row.names = 1)
DEA_6cl_3 <- read.csv("Outputs/Files/Limma/Cluster_3_DEGs_all.csv", row.names = 1)
DEA_6cl_4 <- read.csv("Outputs/Files/Limma/Cluster_4_DEGs_all.csv", row.names = 1)
DEA_6cl_5 <- read.csv("Outputs/Files/Limma/Cluster_5_DEGs_all.csv", row.names = 1)
DEA_6cl_6 <- read.csv("Outputs/Files/Limma/Cluster_6_DEGs_all.csv", row.names = 1)

# 2. Function to filter and add cluster label
filter_deg <- function(df, cluster_name) {
  df %>%
    mutate(gene = rownames(df)) %>%
    filter(adj.P.Val < 0.01, logFC > 0
           ) %>%   # significance filter
    arrange(adj.P.Val) %>%                   # sort most significant → least
   slice_head(n = 100) %>%                  # keep only top 100
    mutate(cluster = cluster_name)
}

# 3. Apply to all clusters
deg1 <- filter_deg(DEA_6cl_1, "Cluster_1")
deg2 <- filter_deg(DEA_6cl_2, "Cluster_2")
deg3 <- filter_deg(DEA_6cl_3, "Cluster_3")
deg4 <- filter_deg(DEA_6cl_4, "Cluster_4")
deg5 <- filter_deg(DEA_6cl_5, "Cluster_5")
deg6 <- filter_deg(DEA_6cl_6, "Cluster_6")

# 4. Combine into one dataframe
degs_6cl <- bind_rows(deg1, deg2, deg3, deg4, deg5, deg6)

write.csv(degs_6cl, file = file.path(output_dir, "DEGs_6cl.csv"), row.names = TRUE)