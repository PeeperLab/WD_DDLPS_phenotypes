# This script performs consensus clustering on the proteomic dataset (WD / DD samples only)

library(openxlsx)
library(dplyr)
library(ConsensusClusterPlus)

set.seed(1507)
## ---- Working directory ----------------------------------------------------
## Before running, set the working directory to the "Proteomics" folder of
## this repository (the folder that holds Inputs/ and Outputs/), e.g.
## setwd("/path/to/WD_DDLPS_phenotypes/Proteomics")

proteomic_data <- read.csv("Inputs/Proteomic_data_batch_corrected.csv", row.names = 1)
proteomic_metadata <- read.xlsx("Inputs/Proteomics_metadata.xlsx")

output_dir <- "./Outputs/Files/Consensus_clustering"
if (!dir.exists(output_dir)) {dir.create(output_dir, recursive = TRUE)}

# Put ID as rownames and order
rownames(proteomic_metadata) <- proteomic_metadata[[1]]   # use first column as rownames

#select tumor samples (exclude healthy fat)
proteomic_metadata <- proteomic_metadata[
  proteomic_metadata$class == "DD" | proteomic_metadata$class == "WD",]

proteomic_data <- proteomic_data[rownames(proteomic_data) %in% rownames(proteomic_metadata), ]

proteomic_data <- t(proteomic_data) #samples should be columns

#If you want to select proteins with most variance or expression:
protein_variance <- apply(proteomic_data, 1, var)
selected_proteins <- names(sort(protein_variance, decreasing = TRUE))[1:500]

proteomic_data_filtered <- proteomic_data[selected_proteins, ]
#proteomic_data_filtered <- scale(proteomic_data_filtered)

results <- ConsensusClusterPlus(
  proteomic_data_filtered,              # Input data matrix
  maxK = 8,                 # Maximum number of clusters to test
  reps = 1000,              # Number of resampling iterations
  pItem = 0.8,              # Proportion of samples for resampling (e.g., 80%)
  pFeature = 1,             # Proportion of genes used (1 = all genes)
  clusterAlg = "km",       # Clustering algorithm (options: "km", "pam", "hc", etc.)
  distance = "euclidean",   # Distance metric (e.g., "pearson", "euclidean")
  seed = 1507,              # Set seed for reproducibility
  plot = "pdf",              # Output format 
  title = output_dir     
)

icluster <- 6  # Select optimal k 
#heatmap(results[[icluster]]$consensusMatrix)
consensus_clusters <- results[[icluster]]$consensusClass

consensus_clusters_df <- as.data.frame(consensus_clusters)

# Re-order clusters to custom order
cluster_map <- c(
  "1" = 1,
  "2" = 3,
  "3" = 4,
  "4" = 2,
  "5" = 6,
  "6" = 5
)

# Create a new data frame with switched cluster labels
consensus_clusters_df_ordered <- consensus_clusters_df

consensus_clusters_df_ordered$consensus_clusters <-
  cluster_map[as.character(consensus_clusters_df_ordered$consensus_clusters)]

write.csv(consensus_clusters_df_ordered, file = file.path(output_dir, "Consensus_clustering.csv"))
