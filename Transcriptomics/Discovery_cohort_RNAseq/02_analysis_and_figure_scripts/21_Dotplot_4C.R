# This script creates a dotplot for the expression of TGFb and MSC related genes across the 6 clusters identified in the consensus clustering analysis.

library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(readr)
library(stringr)  # str_extract()

# set working directory
set.seed(1507)
## ---- Working directory ----------------------------------------------------
## Before running, set the working directory to the "Transcriptomics/Discovery_cohort_RNAseq" folder of
## this repository (the folder that holds Inputs/ and Outputs/), e.g.
## setwd("/path/to/WD_DDLPS_phenotypes/Transcriptomics/Discovery_cohort_RNAseq")

output_dir <- "Outputs/Figures/Dotplot_fig_4C"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

metadata <- read_csv("Inputs/Metadata_discovery_WD_DD.csv")
exprs_data    <- readRDS("Inputs/Batch_corrected_data_discovery_WD_DD.rds")
clusters6 <- read.csv("Outputs/Files/Consensus_clustering/Consensus_clustering.csv", row.names = 1)

metadata <- metadata %>%
  left_join(clusters6, by = "sampleID") %>%
  rename(Cluster = consensus_clusters)

# Subset target clusters (optional)
target_clusters <- c(1,2,3,4,5,6)
metadata_sub <- metadata %>% filter(Cluster %in% target_clusters)
samples_sub <- metadata_sub$sampleID
exprs_sub <- exprs_data[, samples_sub]

# =============================
# Scale expression across samples (z-score)
# =============================
scaled_expr <- t(scale(t(exprs_sub)))  # genes x samples
scaled_expr_df <- as.data.frame(t(scaled_expr)) %>% rownames_to_column("sampleID")
scaled_cluster_df <- left_join(
  scaled_expr_df,
  metadata_sub %>% select(sampleID, Cluster),
  by = "sampleID"
)

# Define target genes
#TGFb & MSC related
target_genes <- c(
  "TGFB1", "TGFBR1", 
  "TGFBI",  "FN1", "SERPINE1", "TGFBR3",
  "TGFB2", "TGFB3",
 "FSTL1", "SMURF2", "FBN1", "PDGFRA", "CD34", "NT5E", "DPP4", "DLK1"
)

scaled_cluster_df <- scaled_cluster_df %>% select(all_of(target_genes), Cluster)

# Pivot to long format
expr_data_long <- scaled_cluster_df %>%
  pivot_longer(cols = all_of(target_genes), names_to = "Gene", values_to = "Expression") %>%
  mutate(
    Gene = factor(Gene, levels = rev(target_genes)),  # preserves top-to-bottom
    Cluster = factor(Cluster, levels = as.character(target_clusters))  # explicit cluster order
  )

# Compute cluster-average expression
summary_df <- expr_data_long %>%
  group_by(Cluster, Gene) %>%
  summarise(avg_expr = mean(Expression, na.rm = TRUE), .groups = "drop")

# Load DEA results for clusters
dea_dir <- "Outputs/Files/Limma"
dea_files <- list.files(dea_dir, pattern = "Cluster_\\d+_DEGs_all.csv", full.names = TRUE)

dea_list <- lapply(dea_files, function(f) {
  df <- read.csv(f, row.names = 1, check.names = FALSE)
  df <- df %>% rownames_to_column("Gene")
  df <- df %>% mutate(
    signif = case_when(
      adj.P.Val < 0.0001 ~ "****",
      adj.P.Val < 0.001  ~ "***",
      adj.P.Val < 0.01   ~ "**",
      adj.P.Val < 0.05   ~ "*",
      TRUE               ~ ""
    )
  )
  return(df)
})

names(dea_list) <- paste0("Cluster", seq_along(dea_list))

# Combine DEA into one table: Cluster | Gene | adj.P.Val | signif
dea_tidy <- bind_rows(
  lapply(names(dea_list), function(cl) {
    df <- dea_list[[cl]]
    df %>% mutate(Cluster = str_extract(cl, "\\d+"))
  })
)

# Keep only target genes
dea_tidy <- dea_tidy %>% filter(Gene %in% target_genes)

# Merge DEA with expression summary
merged_df <- summary_df %>%
  left_join(
    dea_tidy %>% select(Gene, Cluster, adj.P.Val, signif),
    by = c("Gene", "Cluster")
  ) %>%
  mutate(
    neg_log10_adjP = -log10(adj.P.Val),
    neg_log10_adjP = ifelse(is.infinite(neg_log10_adjP), NA, neg_log10_adjP)
  )

# Dotplot
dotplot <- ggplot(
  merged_df,
  aes(x = Cluster, y = Gene, color = avg_expr, size = neg_log10_adjP)
) +
  geom_point() +
  scale_color_gradient2(low = "#80B1D2", mid = "white", high = "darkred", midpoint = 0) +
  scale_size_continuous(name = "-log10(adj.P.Val)", breaks = scales::pretty_breaks(n = 5)) +
  scale_y_discrete(limits = rev(target_genes)) +  # enforce gene order
  scale_x_discrete(limits = as.character(target_clusters)) +  # enforce cluster order
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "RNA-seq: Expression and DEA significance",
    x = "Cluster",
    y = "Gene",
    color = "Mean Expression"
  )

dotplot

# Save plot
ggsave(
  filename = file.path(output_dir, "Dotplot_4C.pdf"),
  plot = dotplot,
  width = 4,
  height = 5,
)