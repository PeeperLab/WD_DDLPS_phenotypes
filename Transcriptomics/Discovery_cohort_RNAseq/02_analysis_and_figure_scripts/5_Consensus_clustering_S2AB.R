### This script is used to perform consensus clustering on the WD/DD subset of the transcriptomic discovery cohort

library(ConsensusClusterPlus)
library(readr)
library(dplyr)
library(matrixStats)
library(tibble)

# set working directory
set.seed(1507)
## ---- Working directory ----------------------------------------------------
## Before running, set the working directory to the "Transcriptomics/Discovery_cohort_RNAseq" folder of
## this repository (the folder that holds Inputs/ and Outputs/), e.g.
## setwd("/path/to/WD_DDLPS_phenotypes/Transcriptomics/Discovery_cohort_RNAseq")

output_dir <- "Outputs/Files/Consensus_clustering"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# load data
exprs_data <- readRDS("Inputs/Batch_corrected_data_discovery_WD_DD.rds")


#selecte the top genes with highest variance
gene_variances <- rowVars(as.matrix(exprs_data))
top_genes <- names(sort(gene_variances, decreasing = TRUE))[1:5000]
exprs_data <- exprs_data[top_genes, ]

count_matrix <- as.matrix(exprs_data)

# ConsensusClusterPlus
results <- ConsensusClusterPlus(
  count_matrix,                        
  maxK = 10,                        
  reps = 1000,                    
  pItem = 0.8,                        
  pFeature = 1,                       
  clusterAlg = "km",                  
  distance = "euclidean",             
  seed = 1507,                        
  plot = "pdf",                       
  title = output_dir                
)

#Save clustering k6
bestK_clusters <- as.data.frame(results[[6]]$consensusClass)
colnames(bestK_clusters) <- "consensus_clusters"
bestK_clusters <- bestK_clusters %>% rownames_to_column(var = "sampleID")

bestK_clusters <- bestK_clusters %>%
  mutate(
    consensus_clusters = case_when(
      consensus_clusters == 2 ~ 3,
      consensus_clusters == 3 ~ 2,
      TRUE ~ consensus_clusters
    )
  )

write.csv(bestK_clusters, file = file.path(output_dir, "Consensus_clustering.csv"))
