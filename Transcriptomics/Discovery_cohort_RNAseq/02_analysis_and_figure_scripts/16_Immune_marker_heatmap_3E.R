# This script creates a heatmap of immune cell marker expression averaged by cluster for clusters 1 and 2, with significance stars from DEA results.

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)
library(RColorBrewer)

# set working directory
set.seed(1507)
## ---- Working directory ----------------------------------------------------
## Before running, set the working directory to the "Transcriptomics/Discovery_cohort_RNAseq" folder of
## this repository (the folder that holds Inputs/ and Outputs/), e.g.
## setwd("/path/to/WD_DDLPS_phenotypes/Transcriptomics/Discovery_cohort_RNAseq")

output_dir <- "Outputs/Figures/Immune_markers_3E"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

metadata <- read_csv("Inputs/Metadata_discovery_WD_DD.csv")
exprs_data    <- readRDS("Inputs/Batch_corrected_data_discovery_WD_DD.rds")
clusters6 <- read.csv("Outputs/Files/Consensus_clustering/Consensus_clustering.csv", row.names = 1)

# Merge clusters into metadata
metadata <- metadata %>%
  left_join(clusters6, by = "sampleID")

colnames(metadata)[ncol(metadata)] <- "clusters6"

# Define immune cell markers
immune_markers <- list(
  "T cells" = c("CD3D", "CD8A", "CD4", "FOXP3", "IL2RA", "GZMA", "GZMB"),
  "B cells" = c("CD19", "CD79A", "CD79B", "MZB1"),
  "Dendritic cells" = c("CLEC9A", "CD1C", "FCER1A", "CLEC10A", "FCER2"),
  "Macrophages (pan)" = c("CD68", "CSF1R", "CSF1","IL1B", "CD80", "CD86", "CD163")
)

# Flatten list for annotation
marker_df <- bind_rows(
  lapply(names(immune_markers), function(cell_type) {
    tibble(Gene = immune_markers[[cell_type]], CellType = cell_type)
  })
)

# Keep only genes present in expression data
marker_df <- marker_df %>% filter(Gene %in% rownames(exprs_data))
target_genes <- marker_df$Gene

sub_exprs_data <- exprs_data[target_genes, ]

# Scale across all samples
scaled_exprs <- t(scale(t(sub_exprs_data)))
scaled_expr_df <- as.data.frame(t(scaled_exprs)) %>% rownames_to_column("sampleID")

# Subset metadata for clusters 1 and 2
sub_metadata <- metadata %>%
  select(sampleID, clusters6) %>%
  filter(clusters6 %in% c("1", "2"))

scaled_cluster_df <- left_join(
  scaled_expr_df,
  sub_metadata,
  by = "sampleID"
) %>%
  filter(!is.na(clusters6))

# Pivot to long format
expr_data_long <- scaled_cluster_df %>%
  select(all_of(target_genes), clusters6) %>%
  pivot_longer(
    cols = all_of(target_genes),
    names_to = "Gene",
    values_to = "Expression"
  ) %>%
  left_join(marker_df, by = "Gene")

# Preserve gene order according to immune_markers
ordered_genes <- unlist(immune_markers, use.names = FALSE)
ordered_genes <- intersect(ordered_genes, rownames(exprs_data))

# Set factor levels in expr_data_long
expr_data_long <- expr_data_long %>%
  mutate(
    Gene = factor(Gene, levels = rev(ordered_genes)),  # reverse for top-to-bottom plotting
    clusters6 = factor(clusters6, levels = c("1", "2"))
  )

# Compute cluster-average expression
summary_df <- expr_data_long %>%
  group_by(clusters6, Gene, CellType) %>%
  summarise(avg_expr = mean(Expression, na.rm = TRUE), .groups = "drop")

# Now fix the gene order
summary_df <- summary_df %>%
  mutate(
    Gene = factor(Gene, levels = rev(ordered_genes)),  # top-to-bottom order
    clusters6 = factor(clusters6, levels = c("1", "2"))
  )

# Define colors for cell types
cell_types <- unique(summary_df$CellType)
cell_colors <- setNames(
  colorRampPalette(brewer.pal(8, "Set2"))(length(cell_types)),
  cell_types
)

# Load DEA results and add significance stars
dea_outdir <- "./Outputs/Files/Limma/Cluster_1_vs_2"
dir.create(dirname(dea_outdir), recursive = TRUE, showWarnings = FALSE)
dea_file <- file.path(dea_outdir, "Cluster_1_vs_2_DEGs.csv")

dea_df <- read.csv(dea_file, row.names = 1, check.names = FALSE) %>%
  rownames_to_column("Gene") %>%
  mutate(
    signif = case_when(
      adj.P.Val < 0.0001 ~ "****",
      adj.P.Val < 0.001  ~ "***",
      adj.P.Val < 0.01   ~ "**",
      adj.P.Val < 0.05   ~ "*",
      TRUE               ~ "ns"
    )
  )

summary_df <- summary_df %>%
  left_join(dea_df %>% select(Gene, adj.P.Val, signif), by = "Gene")

# Right-side position for stars
star_x <- length(levels(summary_df$clusters6)) + 0.7

# Plot heatmap
heatmap <- ggplot(summary_df, aes(x = clusters6, y = Gene, fill = avg_expr)) +
  geom_tile(color = "white") +
  scale_y_discrete(limits = rev(ordered_genes)) +
  scale_fill_gradient2(
    low = "#80B1D2",
    mid = "white",
    high = "darkred",
    midpoint = 0
  ) +
  # Add stars on right
  geom_text(
    data = summary_df %>% group_by(Gene) %>% slice(1) %>% ungroup(),
    aes(y = Gene, x = star_x, label = signif),
    inherit.aes = FALSE,
    size = 4
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    axis.text.y = element_text(color = cell_colors[summary_df$CellType]),
    plot.margin = margin(5.5, 50, 5.5, 5.5)  # extra space for stars
  ) +
  labs(
    title = "Cluster 1 vs 2: Immune cell marker expression",
    subtitle = "Color = mean Z-score (scaled across all samples); stars = DEA significance",
    x = "Cluster",
    y = "Gene",
    fill = "Mean Z-score"
  ) +
  coord_cartesian(clip = "off")

# Show heatmap
heatmap

ggsave(
  filename = file.path(output_dir, "Heatmap_immune_markers.pdf"),
  plot = heatmap,
  width = 5,
  height = 6,
  units = "in")