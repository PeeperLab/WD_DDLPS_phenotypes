# This script creates dot plots of selected proteins per cluster in the proteomic cohort.
# Colour: mean scaled protein expression per cluster; size: -log10 adjusted p-value of the
# one-versus-rest limma analysis (8_Limma_clusters_tumor.R).
# Fig. 3F: immune-related proteins; Fig. 4D: TGF-beta / MSC-related proteins.

library(openxlsx)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(ggplot2)
library(scales)

set.seed(1507)
## ---- Working directory ----------------------------------------------------
## Before running, set the working directory to the "Proteomics" folder of
## this repository (the folder that holds Inputs/ and Outputs/), e.g.
## setwd("/path/to/WD_DDLPS_phenotypes/Proteomics")

output_dir <- "Outputs/Figures/Dotplot_markers_3F_4D"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Load data
proteomic_data <- read.csv("Inputs/Proteomic_data_batch_corrected.csv", row.names = 1)
proteomic_metadata <- read.xlsx("Inputs/Proteomics_metadata.xlsx")
rownames(proteomic_metadata) <- proteomic_metadata[[1]]

# Tumour samples only (DD and WD), as used for consensus clustering
proteomic_metadata <- proteomic_metadata[
  proteomic_metadata$class == "DD" | proteomic_metadata$class == "WD", ]
proteomic_data <- proteomic_data[rownames(proteomic_data) %in% rownames(proteomic_metadata), ]

# Add consensus clusters, matched by sample ID
clusters <- read.csv("Outputs/Files/Consensus_clustering/Consensus_clustering.csv", row.names = 1)
proteomic_metadata$Cluster <- clusters$consensus_clusters[match(rownames(proteomic_metadata), rownames(clusters))]
stopifnot(!anyNA(proteomic_metadata$Cluster))

# Scale the expression data (z-score per protein across samples)
scaled_expr <- as.data.frame(scale(proteomic_data))
scaled_expr$ID <- rownames(proteomic_data)
scaled_cluster_df <- left_join(scaled_expr, proteomic_metadata, by = "ID")

cluster_order <- c("1", "2", "3", "4", "5", "6")
scaled_cluster_df$Cluster <- factor(scaled_cluster_df$Cluster, levels = cluster_order)

# limma results (one file per cluster): Cluster | Protein | adj.P.Val
files <- list.files("Outputs/Files/Limma/", pattern = "*.csv", full.names = TRUE)
limma_list <- lapply(files, read.csv, row.names = 1)
names(limma_list) <- gsub("Outputs/Files/Limma/|_DE_proteins.csv", "", files)
limma_tidy <- bind_rows(
  lapply(names(limma_list), function(nm) {
    limma_list[[nm]] %>%
      rownames_to_column("Protein") %>%
      dplyr::select(Protein, adj.P.Val) %>%
      mutate(Cluster = str_extract(nm, "[0-9]+"))
  })
)

# Proteins shown in each figure panel
panels <- list(
  Fig3F_immune = c(
    "PTPRC", "CD48", "MZB1", "CD38", "GZMK", "CD74",
    "ITGAM", "MMP9", "LYZ", "ARG1", "S100A8", "S100A9", "ELANE", "MPO",
    "FCN1", "POSTN", "SPP1", "APOE",
    "PLAUR", "ITGAX",
    "TYROBP",
    "STAT1", "PSME3",
    "GBP5",
    "MSR1", "TAP1", "LGALS3", "CTSS", "CTSD", "VSIG4",
    "FCGR1A", "FCGR3A",
    "CSF1R", "CD163",
    "AIF1",
    "MRC1",
    "CD209",
    "FOLR2", "IL18", "CXCL12",
    "CPA3", "CMA1", "CTSG"
  ),
  Fig4D_TGFb_MSC = c(
    "TGFB1", "TGFBI", "FN1", "SERPINE1",
    "TGFBR3", "FSTL1", "SMURF2", "PDGFRA", "CD34", "NT5E", "DPP4"
  )
)
plot_heights <- c(Fig3F_immune = 6.5, Fig4D_TGFb_MSC = 3.5)

make_dotplot <- function(target_genes) {
  # Mean scaled expression per cluster and protein
  summary_df <- scaled_cluster_df %>%
    dplyr::select(all_of(target_genes), Cluster) %>%
    pivot_longer(cols = all_of(target_genes), names_to = "Protein", values_to = "Expression") %>%
    group_by(Cluster, Protein) %>%
    summarise(avg_expr = mean(Expression, na.rm = TRUE), .groups = "drop")

  # Add limma significance
  merged_df <- summary_df %>%
    left_join(limma_tidy %>% filter(Protein %in% target_genes), by = c("Cluster", "Protein")) %>%
    mutate(
      neg_log10_adjP = -log10(adj.P.Val),
      neg_log10_adjP = ifelse(is.infinite(neg_log10_adjP), NA, neg_log10_adjP)
    )
  merged_df$Protein <- factor(merged_df$Protein, levels = rev(target_genes))
  merged_df$Cluster <- factor(merged_df$Cluster, levels = cluster_order)

  ggplot(merged_df, aes(x = Cluster, y = Protein, color = avg_expr, size = neg_log10_adjP)) +
    geom_point() +
    scale_color_gradient2(low = "#80B1D2", mid = "white", high = "darkred", midpoint = 0) +
    scale_size_continuous(name = "-log10(adj.P.Val)", breaks = scales::pretty_breaks(n = 5)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      title = "Dotplot: Expression (color) and LIMMA significance (size)",
      x = "Cluster",
      y = "Protein",
      color = "Mean Expression"
    )
}

for (panel in names(panels)) {
  dotplot <- make_dotplot(panels[[panel]])
  print(dotplot)
  ggsave(
    filename = file.path(output_dir, paste0("Dotplot_", panel, ".pdf")),
    plot = dotplot,
    width = 4,
    height = plot_heights[[panel]],
    dpi = 300
  )
}
