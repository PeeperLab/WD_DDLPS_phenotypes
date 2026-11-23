# This script creates a UMAP of the WD/DD subset of the transcriptomic discovery cohort colored by the consensus clusters obtained from consensus clustering.

library(dplyr)
library(ggplot2)
library(readr)
library(ggrepel)
library(umap)

# set working directory
set.seed(1507)
## ---- Working directory ----------------------------------------------------
## Before running, set the working directory to the "Transcriptomics/Discovery_cohort_RNAseq" folder of
## this repository (the folder that holds Inputs/ and Outputs/), e.g.
## setwd("/path/to/WD_DDLPS_phenotypes/Transcriptomics/Discovery_cohort_RNAseq")

output_dir <- "Outputs/Figures/UMAP_clusters_2A"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# load data
metadata <- read_csv("Inputs/Metadata_discovery_WD_DD.csv")
data <- readRDS("Inputs/Batch_corrected_data_discovery_WD_DD.rds")
clusters <- read.csv("Outputs/Files/Consensus_clustering/Consensus_clustering.csv", row.names = 1)

gene_vars <- apply(data, 1, var)
top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:5000]
data_hvg <- data[top_genes, ]

# set UMAP parameters
umap_config <- umap.defaults
umap_config$n_neighbors <- 15 
umap_config$min_dist <- 0.9  

# UMAP
umap_result <- umap(t(data_hvg), config = umap_config)

umap_data <- as.data.frame(umap_result$layout)
colnames(umap_data) <- c("UMAP1", "UMAP2")

# match metadata rows to the samples (columns) of the expression matrix by sample ID
metadata <- metadata[match(colnames(data_hvg), metadata$sampleID), ]
stopifnot(identical(metadata$sampleID, colnames(data_hvg)))
umap_data$subtype <- metadata$subtype
umap_data$dataset <- metadata$dataset
umap_data$sampleID <- metadata$sampleID
umap_data$dataset <- as.factor(umap_data$dataset)

clusters$sampleID <- rownames(clusters)
umap_data$sampleID <- rownames(umap_data)
umap_data$consensus_clusters <- clusters$consensus_clusters

umap_data$consensus_clusters <- as.factor(umap_data$consensus_clusters)

cluster_colors1 <- c(
  "1" = "#A5C1A4",
  "2" = "#827BB3",
  "3" = "#DF9677",
  "4" = "#1F70B5",
  "5" = "#F7CBDF",
  "6" = "#FDEDB9"
)

umap2 <- ggplot(umap_data, aes(x = UMAP1, y = UMAP2, color = consensus_clusters, label = NA)) +
  geom_point(size = 5.5, alpha = 0.8, shape = 16) +
  coord_fixed() +
  theme_minimal() +
  labs(title = "", color = "Sample type") +
  theme(legend.position = "right",
        aspect.ratio = 1,
        panel.background = element_rect(fill = "white",colour = NA),
        plot.background = element_rect(fill = "white",colour = NA),
        text = element_text(size = 16), 
        axis.text = element_text(size = 16), 
        axis.title = element_text(size = 16, face = "bold"),  
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5), 
        legend.text = element_text(size = 16),  
        legend.title = element_text(size = 16),
        panel.grid = element_blank())+
  scale_color_manual(
    values = cluster_colors1,
    name = "Cluster"
  ) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

umap2

ggsave(
  filename = file.path(output_dir, "RNAseq_umap_clusters.pdf"),
  plot = umap2,
  width = 7.5,
  height = 6,
  units = "in")