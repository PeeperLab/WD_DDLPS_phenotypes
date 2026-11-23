# This script creates a volcano plot for each cluster based on DEA between all clusters

library(tidyverse)
library(ggplot2)
library(dplyr)
library(readr)
library(ggrepel)

# set working directory
set.seed(1507)
## ---- Working directory ----------------------------------------------------
## Before running, set the working directory to the "Transcriptomics/Discovery_cohort_RNAseq" folder of
## this repository (the folder that holds Inputs/ and Outputs/), e.g.
## setwd("/path/to/WD_DDLPS_phenotypes/Transcriptomics/Discovery_cohort_RNAseq")

output_dir <- "Outputs/Figures/Volcano_plots_fig_4A_5A"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Load cluster DEA results
cluster_files <- list.files(
  "./Outputs/Files/Limma",
  pattern = "Cluster_\\d+_DEGs_all.csv",
  full.names = TRUE
)

dea_results_list <- lapply(cluster_files, read.csv, row.names = 1)
names(dea_results_list) <- paste0("Cluster", seq_along(dea_results_list))

# Volcano plot parameters
logFC_cutoff <- 1
padj_cutoff  <- 0.05

# Loop over clusters
for (cluster_name in names(dea_results_list)) {
  
  cat("Making volcano plot for", cluster_name, "...\n")
  
  res <- dea_results_list[[cluster_name]]
  
  # Prepare dataframe
  volcano_df <- res %>%
    rownames_to_column(var = "gene") %>%
    mutate(
      neg_log10_padj = -log10(adj.P.Val),
      significance = case_when(
        adj.P.Val < padj_cutoff & logFC >= logFC_cutoff  ~ "Upregulated",
        adj.P.Val < padj_cutoff & logFC <= -logFC_cutoff ~ "Downregulated",
        TRUE                                             ~ "Not significant"
      )
    )
  
  top_genes <- volcano_df %>%
    filter(adj.P.Val < padj_cutoff) %>%
    arrange(adj.P.Val) %>%
    slice_head(n = 20)
  

  # Volcano plot
  p_volcano <- ggplot(
    volcano_df,
    aes(x = logFC, y = neg_log10_padj)
  ) +
    geom_point(
      aes(color = significance),
      alpha = 0.7,
      size = 1.5,
      shape = 16
    ) +
    geom_text_repel(
      data = top_genes,
      aes(label = gene),
      size = 3.5,
      max.overlaps = Inf,
      box.padding = 0.4,        
      point.padding = 0.3,     
      segment.color = "grey50",
      segment.size = 0.3
    )+
    scale_color_manual(
      values = c(
        "Upregulated"   = "darkred",
        "Downregulated" = "#80B1D2",
        "Not significant" = "grey70"
      )
    ) +
    geom_vline(
      xintercept = c(-logFC_cutoff, logFC_cutoff),
      linetype = "dashed",
      linewidth = 0.4
    ) +
    geom_hline(
      yintercept = -log10(padj_cutoff),
      linetype = "dashed",
      linewidth = 0.4
    ) +
    theme_minimal() +
    labs(
      title = paste(cluster_name, "Volcano Plot"),
      x = "log2 Fold Change",
      y = expression(-log[10]("Adjusted p-value")),
      color = "Significance"
    ) +
    theme(
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 9, face = "bold"),
      axis.text = element_text(size = 8),
      legend.title = element_text(size = 8, face = "bold"),
      legend.text = element_text(size = 8),
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white"),
      panel.grid = element_blank() 
    )
  
  # Save plot
  ggsave(
    filename = file.path(output_dir, paste0(cluster_name, "_volcano.pdf")),
    plot = p_volcano,
    width = 5,
    height = 4,
    units = "in"
  )

  
  print(p_volcano)
  
  cat("Finished", cluster_name, "\n\n")
}

