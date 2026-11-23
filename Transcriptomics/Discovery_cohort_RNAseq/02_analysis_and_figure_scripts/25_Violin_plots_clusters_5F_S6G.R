# This script creates violin plots of gene expression for a specific gene across the 6 clusters, 
# with significance stars based on DEA results. 

library(readr)
library(dplyr)
library(limma)
library(tibble)
library(stringr)
library(ggplot2)
library(tidyr)

# set working directory
set.seed(1507)
## ---- Working directory ----------------------------------------------------
## Before running, set the working directory to the "Transcriptomics/Discovery_cohort_RNAseq" folder of
## this repository (the folder that holds Inputs/ and Outputs/), e.g.
## setwd("/path/to/WD_DDLPS_phenotypes/Transcriptomics/Discovery_cohort_RNAseq")

output_dir <- "Outputs/Figures/Expression_violins_5F_S6G"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

metadata <- read_csv("Inputs/Metadata_discovery_WD_DD.csv")
exprs_data    <- readRDS("Inputs/Batch_corrected_data_discovery_WD_DD.rds")
clusters6 <- read.csv("Outputs/Files/Consensus_clustering/Consensus_clustering.csv", row.names = 1)

# Merge clusters into metadata
metadata <- metadata %>%
  left_join(clusters6, by = "sampleID")
colnames(metadata)[ncol(metadata)] <- "clusters6"

# Cluster colors
cluster_colors1 <- c(
  "1" = "#A5C1A4",
  "2" = "#827BB3",
  "3" = "#DF9677",
  "4" = "#1F70B5",
  "5" = "#F7CBDF",
  "6" = "#FDEDB9"
)

# Load cluster DEA results
dea_files <- list.files(
  "./Outputs/Files/Limma",
  pattern = "Cluster_\\d+_DEGs_all.csv",
  full.names = TRUE
)

dea_list <- lapply(dea_files, function(f) {
  read.csv(f, row.names = 1, check.names = FALSE) %>%
    rownames_to_column("Gene") %>%
    mutate(
      signif = case_when(
        adj.P.Val < 0.0001 ~ "****",
        adj.P.Val < 0.001  ~ "***",
        adj.P.Val < 0.01   ~ "**",
        adj.P.Val < 0.05   ~ "*",
        TRUE               ~ "ns"
      ),
      Cluster = str_extract(basename(f), "\\d+")
    )
})

dea_all <- bind_rows(dea_list)

# Fig. 5F: ACAN, COL2A1, PROM1; Fig. S6G: SOX2, SOX9, RUNX2
for (gene in c("ACAN", "COL2A1", "PROM1", "SOX2", "SOX9", "RUNX2")) {
  stopifnot(gene %in% rownames(exprs_data))

  # Expression dataframe
  gene_df <- data.frame(
    sampleID   = colnames(exprs_data),
    expression = as.numeric(exprs_data[gene, ])
  ) %>% left_join(metadata, by = "sampleID")

  gene_df$clusters6 <- factor(gene_df$clusters6,
                              levels = sort(unique(gene_df$clusters6)))

  # Filter DEA results for the gene
  dea_gene <- dea_all %>%
    filter(Gene == gene) %>%
    mutate(Cluster = factor(Cluster, levels = names(cluster_colors1)))

  # Position stars slightly above max expression
  max_y <- max(gene_df$expression, na.rm = TRUE)
  dea_gene <- dea_gene %>% mutate(y_pos = max_y * 1.15)

  # Violin plot with statistics
  p <- ggplot(gene_df, aes(x = clusters6, y = expression, fill = clusters6)) +
  
    # Violin
    geom_violin(trim = FALSE, alpha = 0.8) +
  
    # Boxplot
    geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.6) +
  
    # Individual samples
    geom_jitter(width = 0.15, size = 1.6, alpha = 0.7, color = "black", shape = 16) +
  
    # DEA significance stars
    geom_text(
      data = dea_gene,
      aes(x = Cluster, y = y_pos, label = signif),
      inherit.aes = FALSE,
      size = 5
    ) +
  
    # Cluster colors
    scale_fill_manual(values = cluster_colors1) +
  
    theme_minimal() +
    labs(
      title = paste0(gene, " expression per cluster"),
      x = "Cluster",
      y = "Expression"
    ) +
    theme(
      legend.position = "none",
      text = element_text(size = 14),
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      axis.title = element_text(face = "bold")
    )

  print(p)

  ggsave(
    filename = file.path(output_dir, paste0(gene, "_expression.pdf")),
    plot = p,
    width = 6,
    height = 5,
    units = "in")
}
