# This script creates UMAPS of the proteomic dataset (WD / DD samples only) colored by consensus cluster and 
# whether a specific protein is detected in the raw data. 
# It also includes a bar plot quantifying the detection rate of the protein across clusters.

library(readxl)
library(openxlsx)
library(ggplot2)
library(dplyr)
library(tidyr)
library(umap)
library(tibble)  # rownames_to_column(), column_to_rownames()

set.seed(1507)
## ---- Working directory ----------------------------------------------------
## Before running, set the working directory to the "Proteomics" folder of
## this repository (the folder that holds Inputs/ and Outputs/), e.g.
## setwd("/path/to/WD_DDLPS_phenotypes/Proteomics")

proteomic_data <- read.csv("Inputs/Proteomic_data_batch_corrected.csv", row.names = 1)
proteomic_metadata <- read.xlsx("Inputs/Proteomics_metadata.xlsx")
consensus_clusters <- read.csv("Outputs/Files/Consensus_clustering/Consensus_clustering.csv", row.names = 1)

output_dir <- "Outputs/Figures/UMAP_WD_DD_samples_fig_3A"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

output_dir_2 <- "Outputs/Figures/Protein_detection_S6C_F"
if (!dir.exists(output_dir_2)) dir.create(output_dir_2, recursive = TRUE)

rownames(proteomic_metadata) <- proteomic_metadata[[1]] 

#select tumor samples (exclude healthy fat)
proteomic_metadata <- proteomic_metadata[
  proteomic_metadata$class == "DD" | proteomic_metadata$class == "WD",]

proteomic_data <- proteomic_data[rownames(proteomic_data) %in% rownames(proteomic_metadata), ]

# Compute variance and select top 500 proteins
protein_variance <- apply(proteomic_data, 2, var, na.rm = TRUE)  # variance across samples
top_proteins <- names(sort(protein_variance, decreasing = TRUE))[1:500]

# Filter to those top proteins
proteomic_data_filtered <- proteomic_data[, top_proteins]

# Run UMAP
set.seed(1507)
umap_result <- umap(proteomic_data_filtered, n_neighbors = 20, min_dist = 0.8)

# Extract coordinates
umap_coords <- as.data.frame(umap_result$layout)
colnames(umap_coords) <- c("UMAP1", "UMAP2")

# Keep sample IDs
umap_coords$SampleID <- rownames(proteomic_data_filtered)

# Add metadata (consensus cluster, etc.)
umap_coords <- umap_coords %>%
  left_join(proteomic_metadata, by = c("SampleID" = "ID"))

cluster_colors <- c(
  "1" = "#A5C1A4",
  "2" = "#827BB3",
  "3" = "#E34363",
  "4" = "#DF9677",
  "5" = "#1F70B5",
  "6" = "#FDEDB9"
)

# Plot UMAP colored by consensus cluster (Fig 3A)

# Make sure consensus_cluster is a factor
# (matched by sample ID)
umap_coords$consensus_cluster <- as.factor(consensus_clusters[umap_coords$SampleID, "consensus_clusters"])
stopifnot(!anyNA(umap_coords$consensus_cluster))

# Create the UMAP plot
umap_plot <- ggplot(umap_coords, aes(x = UMAP1, y = UMAP2, color = consensus_cluster)) +
  geom_point(size = 5, alpha = 0.8, shape = 16) +
  coord_fixed() +
  scale_color_manual(values = cluster_colors, name = "Consensus Cluster") +
  theme_minimal() +
  labs(
    title = "UMAP of Proteomic Data by Consensus Cluster",
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  theme(
    text = element_text(size = 16),
    aspect.ratio = 1,
    axis.title = element_text(size = 16),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.background = element_rect(fill = "white"),
    panel.grid = element_blank(),
    plot.background = element_rect(fill = "white")
  )

umap_plot

#Save as PDF
ggsave(
  filename = file.path(output_dir, "UMAP_Consensus_Clusters.pdf"),
  plot = umap_plot,
  width = 7,
  height = 5,
  units = "in"
)

# Plot whether protein is detected in raw data (Fig. S6C, S6E)

proteomic_data_raw <- read_excel("Inputs/Raw_combined_data.xlsx")
# Fig. S6C-D (GLI1) and Fig. S6E-F (HHIP)
for (protein_of_interest in c("GLI1", "HHIP")) {

  # Extract detection status for protein of interest from raw data
  detection_df <- proteomic_data_raw %>%
    dplyr::select(ID, all_of(protein_of_interest)) %>%
    setNames(c("SampleID", "raw_value")) %>%
    mutate(Detected = ifelse(is.na(raw_value), "Not detected", "Detected"))

  # Merge with UMAP coordinates
  umap_detection <- umap_coords %>%
    left_join(detection_df, by = "SampleID")

  # Plot
  p_detection <- ggplot(umap_detection, aes(x = UMAP1, y = UMAP2, color = Detected, label = NA)) +
    geom_point(size = 5, alpha = 0.8, shape = 16) +
    geom_text(vjust = -1.2, size = 3) +
    coord_fixed() +
    scale_color_manual(
      values = c("Detected" = "darkred", "Not detected" = "#80B1D2"),
      name = paste0(protein_of_interest, "\nDetection")
    ) +
    theme_minimal() +
    labs(
      x = "UMAP 1",
      y = "UMAP 2"
    ) +
    theme(
      text = element_text(size = 12),
      aspect.ratio = 1,
      axis.title = element_text(size = 12),
      legend.title = element_text(size = 8, face = "bold"),
      legend.text = element_text(size = 8),
      panel.border = element_blank(),
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  p_detection

  ggsave(
    filename = file.path(output_dir_2, paste0("UMAP_", protein_of_interest, "_detection.pdf")),
    plot = p_detection,
    width = 7,
    height = 5,
    units = "in",
  )

  # Quantification of detection rate (bar plots) (Fig. S6D, S6F)

  # Extract detection status for protein of interest from raw data
  detection_df <- proteomic_data_raw %>%
    dplyr::select(ID, all_of(protein_of_interest)) %>%
    setNames(c("SampleID", "raw_value")) %>%
    mutate(Detected = ifelse(is.na(raw_value), "Not detected", "Detected"))

  # Merge with cluster information
  detection_clusters <- consensus_clusters %>%
    rownames_to_column("SampleID") %>%
    left_join(detection_df, by = "SampleID")

  # Compute detection rate per cluster
  detection_rate <- detection_clusters %>%
    group_by(consensus_clusters, Detected) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(consensus_clusters) %>%
    mutate(
      pct = n / sum(n) * 100,
      Detected = factor(Detected, levels = c("Detected", "Not detected"))
    )

  # Create binary grouping: cluster 5 vs all others
  detection_clusters_binary <- detection_clusters %>%
    mutate(cluster_group = ifelse(consensus_clusters == 5, "Cluster 5", "Other clusters"))

  # Create contingency table
  contingency_table_binary <- detection_clusters_binary %>%
    group_by(cluster_group, Detected) %>%
    summarise(n = n(), .groups = "drop") %>%
    pivot_wider(names_from = Detected, values_from = n, values_fill = 0) %>%
    column_to_rownames("cluster_group")

  # Fisher's exact test
  fisher_result <- fisher.test(contingency_table_binary)
  fisher_result

  # Extract p-value
  p_val <- round(fisher_result$p.value, 12)
  p_label <- paste0("p = ", p_val)

  # Plot
  p_detection_bar <- ggplot(detection_rate, aes(x = factor(consensus_clusters), y = pct, fill = Detected)) +
    geom_bar(stat = "identity", position = position_stack(reverse = TRUE)) +
    scale_fill_manual(
      values = c("Detected" = "darkred", "Not detected" = "grey"),
      name = paste0(protein_of_interest, "\nDetection")
    ) +
    annotate("text", x = Inf, y = Inf, label = p_label,
             hjust = 1.1, vjust = 2, size = 3.5) +
    theme_minimal() +
    labs(
      x = "Cluster",
      y = "Percentage of samples (%)",
    ) +
    theme(
      text = element_text(size = 12),
      axis.title = element_text(size = 12),
      legend.title = element_text(size = 8, face = "bold"),
      legend.text = element_text(size = 8),
      panel.border = element_blank(),
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  p_detection_bar

  ggsave(
    filename = file.path(output_dir_2, paste0("Barplot_", protein_of_interest, "_detection.pdf")),
    plot = p_detection_bar,
    width = 7,
    height = 5,
    units = "in",
    device = "pdf"
  )
}
