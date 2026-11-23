# This script creates a heatmap of ssGSEA scores averaged by cluster for the 6 clusters identified in the consensus clustering analysis. 

library(readr)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(pheatmap)
library(grid)
library(tidyr)

# set working directory
set.seed(1507)
## ---- Working directory ----------------------------------------------------
## Before running, set the working directory to the "Transcriptomics/Discovery_cohort_RNAseq" folder of
## this repository (the folder that holds Inputs/ and Outputs/), e.g.
## setwd("/path/to/WD_DDLPS_phenotypes/Transcriptomics/Discovery_cohort_RNAseq")

output_dir <- "Outputs/Figures/ssGSEA_heatmap_2E"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

metadata <- read_csv("Inputs/Metadata_discovery_WD_DD.csv")
k_6clusters <- read.csv("Outputs/Files/Consensus_clustering/Consensus_clustering.csv", row.names = 1)
ssGSEA_res <- read.csv("Outputs/Files/ssGSEA/gsva_results.csv", row.names = 1)

rownames(ssGSEA_res) <- ssGSEA_res$pathway
ssGSEA_res <- ssGSEA_res[,-1]

## add the consensus clustering result to meta data
metadata <- left_join(metadata,k_6clusters, by = "sampleID")
colnames(metadata)[5] <- "clusters6"

# Remove "REACTOME_" prefix from row names
rownames(ssGSEA_res) <- gsub("^REACTOME_", "", rownames(ssGSEA_res))

# Format all other row names (except already customized ones)
format_row_names <- function(name) {
  if (grepl(" ", name)) return(name)  
  name <- gsub("-", " ", name) 
  name <- tolower(name)
  words <- strsplit(name, "_")[[1]]
  words <- sapply(words, function(word) {
    paste0(toupper(substring(word, 1, 1)), substring(word, 2))
  })
  return(paste(words, collapse = " "))
}
rownames(ssGSEA_res) <- sapply(rownames(ssGSEA_res), format_row_names)

# Make sure sample IDs match exactly
intersecting_samples <- intersect(colnames(ssGSEA_res), metadata$sampleID)

# Subset ssGSEA_res to only samples that exist in metadata
ssgsea_sub <- ssGSEA_res[, intersecting_samples, drop = FALSE]

# Match metadata order to ssGSEA columns
metadata_sub <- metadata[match(colnames(ssgsea_sub), metadata$sampleID), ]

# Average ssGSEA by cluster 
ssgsea_long <- as.data.frame(t(ssgsea_sub))
ssgsea_long$sampleID <- rownames(ssgsea_long)

ssgsea_long <- ssgsea_long %>%
  pivot_longer(
    cols = -sampleID,
    names_to = "pathway",
    values_to = "score"
  ) %>%
  left_join(metadata_sub[, c("sampleID", "clusters6")], by = "sampleID")

mean_scores <- ssgsea_long %>%
  group_by(clusters6, pathway) %>%
  summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    names_from = clusters6,
    values_from = mean_score
  )

# Convert back to matrix for heatmap
pathway_matrix <- as.matrix(mean_scores[, -1])
rownames(pathway_matrix) <- mean_scores$pathway

# Keep pathway order
desired_order <- rownames(ssGSEA_res)
pathway_matrix_ordered <- pathway_matrix[desired_order, , drop = FALSE]

# Order clusters 1–6
desired_col_order <- c("1", "2", "3", "4", "5", "6")
pathway_matrix_ordered <- pathway_matrix_ordered[, desired_col_order, drop = FALSE]

annotation_col <- data.frame(
  Cluster = factor(
    colnames(pathway_matrix_ordered),
    levels = c("1", "2", "3", "4", "5", "6")
  )
)

rownames(annotation_col) <- colnames(pathway_matrix_ordered)

annotation_colors <- list(
  Cluster = c(
    "1" = "#A5C1A4",
    "2" = "#827BB3",
    "3" = "#DF9677",
    "4" = "#1F70B5",
    "5" = "#F7CBDF",
    "6" = "#FDEDB9"
  )
)

ht_clusters <- pheatmap(
  pathway_matrix_ordered,
  scale = "row",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  fontsize = 12,
  angle_col = 45,
  color = colorRampPalette(c("#80B1D2", "white", "darkred"))(100),
  border_color = NA,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors
)

ht_clusters

out_file <- file.path(output_dir, "ssGSEA_heatmap.pdf")
# Open PDF device
pdf(out_file, width = 7, height = 6)

print(ht_clusters)

dev.off()
