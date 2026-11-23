# This script creates violin plots of the inferred immune cell fractions from EPIC for the 6 clusters identified in the consensus clustering analysis, 
# and a heatmap of logFC with significance stars for pairwise cluster comparisons.

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(EPIC)
library(rstatix)
library(scales)
library(tibble)

# set working directory
set.seed(1507)
## ---- Working directory ----------------------------------------------------
## Before running, set the working directory to the "Transcriptomics/Discovery_cohort_RNAseq" folder of
## this repository (the folder that holds Inputs/ and Outputs/), e.g.
## setwd("/path/to/WD_DDLPS_phenotypes/Transcriptomics/Discovery_cohort_RNAseq")

output_dir <- "Outputs/Figures/Immune_infiltration_3D_S5C"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

metadata <- read_csv("Inputs/Metadata_discovery_WD_DD.csv")
exprs_data    <- readRDS("Inputs/Batch_corrected_data_discovery_WD_DD.rds")
clusters6 <- read.csv("Outputs/Files/Consensus_clustering/Consensus_clustering.csv", row.names = 1)

# Merge cluster info
metadata <- metadata %>%
  left_join(clusters6, by = "sampleID") %>%
  rename(clusters6 = consensus_clusters)

# Run EPIC for all clusters
epic_res <- EPIC(bulk = exprs_data)

# Extract cell-type fractions
fractions <- epic_res$cellFractions  # rows = samples, cols = cell types

# Convert to long format for plotting
fractions_long <- as.data.frame(fractions) %>%
  rownames_to_column("sampleID") %>%
  left_join(metadata %>% select(sampleID, clusters6), by = "sampleID") %>%
  pivot_longer(
    cols = -c(sampleID, clusters6),
    names_to = "CellType",
    values_to = "Fraction"
  )

# Average fractions per cluster
summary_df <- fractions_long %>%
  group_by(clusters6, CellType) %>%
  summarise(avg_fraction = mean(Fraction, na.rm = TRUE), .groups = "drop")

# Plot violin with dynamic y-axis and p-value stars
violins <- ggplot(
  fractions_long,
  aes(x = factor(clusters6), y = Fraction, fill = factor(clusters6))
) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_jitter(
    width = 0.15,
    size = 0.8,
    alpha = 0.7,
    color = "black",
    shape = 16
  ) +
  facet_wrap(~ CellType, scales = "free_y") +
  scale_fill_manual(
    values = c(
      "1" = "#A5C1A4",
      "2" = "#827BB3",
      "3" = "#DF9677",
      "4" = "#1F70B5",
      "5" = "#F7CBDF",
      "6" = "#FDEDB9"
    )
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    strip.text = element_text(face = "bold"),
    panel.border = element_rect(color = "black", fill = NA)
  ) +
  labs(
    x = "Cluster",
    y = "EPIC estimated fraction"
  )

violins

ggsave(
  filename = file.path(output_dir, "Violin_plots.pdf"),
  plot = violins,
  width = 8,
  height = 6,
  units = "in")

# Statistics
# Ensure clusters are factors
fractions_long <- fractions_long %>%
  mutate(clusters6 = factor(clusters6, levels = 1:6))

# Compute mean fraction per cluster and cell type
mean_df <- fractions_long %>%
  group_by(CellType, clusters6) %>%
  summarise(mean_frac = mean(Fraction, na.rm = TRUE), .groups = "drop")

# Generate lower triangle pairs: cluster1 < cluster2
cluster_levels <- levels(fractions_long$clusters6)
lower_pairs <- expand.grid(cluster1 = cluster_levels,
                           cluster2 = cluster_levels,
                           stringsAsFactors = FALSE) %>%
  filter(as.numeric(cluster1) < as.numeric(cluster2))

# Compute logFC
logfc_df <- lower_pairs %>%
  left_join(mean_df %>% rename(cluster1 = clusters6, mean1 = mean_frac),
            by = "cluster1") %>%
  left_join(mean_df %>% rename(cluster2 = clusters6, mean2 = mean_frac),
            by = c("cluster2", "CellType")) %>%
  mutate(logFC = log2(mean1 + 1e-6) - log2(mean2 + 1e-6))

# Perform Dunn tests for the lower triangle
dunn_results <- fractions_long %>%
  group_by(CellType) %>%
  dunn_test(Fraction ~ clusters6, p.adjust.method = "BH") %>%
  filter(as.numeric(group1) < as.numeric(group2)) %>%
  rename(cluster1 = group1, cluster2 = group2)

# Combine logFC and statistics
heatmap_df_tri <- logfc_df %>%
  left_join(dunn_results, by = c("CellType", "cluster1", "cluster2")) %>%
  mutate(signif = case_when(
    p.adj < 0.0001 ~ "****",
    p.adj < 0.001  ~ "***",
    p.adj < 0.01   ~ "**",
    p.adj < 0.05   ~ "*",
    TRUE           ~ ""
  ))


# Bottom-right triangle:
# Flip only x-axis
heatmap_df_tri_plot <- heatmap_df_tri %>%
  mutate(
    x_plot = factor(cluster2, levels = 6:1),  # reverse x-axis
    y_plot = factor(cluster1, levels = 1:6)   # keep y-axis normal
  )


# Rescale logFC per CellType
heatmap_df_rescaled <- heatmap_df_tri_plot %>%
  group_by(CellType) %>%
  mutate(
    logFC_scaled = scales::rescale(logFC, to = c(-1, 1))  # scale to [-1, 1] for gradient
  ) %>%
  ungroup()

# Plot
heatmap_stats <- ggplot(heatmap_df_rescaled, aes(x = x_plot, y = y_plot)) +
  geom_tile(aes(fill = logFC_scaled), color = "white") +
  geom_text(aes(label = signif), size = 4) +
  facet_wrap(~ CellType, ncol = 4, scales = "free") +
  scale_fill_gradient2(
    low = "#80B1D2",
    mid = "white",
    high = "darkred",
    midpoint = 0,
    name = "log2 FC (scaled)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank()
  ) +
  labs(
    title = "EPIC immune cell fractions",
    subtitle = "Bottom-right triangle: log2FC with BH-adjusted significance\nEach cell type scaled independently"
  )

heatmap_stats

ggsave(
  filename = file.path(output_dir, "Heatmap_stats.pdf"),
  plot = heatmap_stats,
  width = 9,
  height = 6,
  units = "in")