# This script creates a heatmap of the top 100 upregulated genes in each of the 6 clusters obtained from 
# consensus clustering on the WD/DD subset of the transcriptomic discovery cohort.

library(readr)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(grid)

# set working directory
set.seed(1507)
## ---- Working directory ----------------------------------------------------
## Before running, set the working directory to the "Transcriptomics/Discovery_cohort_RNAseq" folder of
## this repository (the folder that holds Inputs/ and Outputs/), e.g.
## setwd("/path/to/WD_DDLPS_phenotypes/Transcriptomics/Discovery_cohort_RNAseq")

output_dir <- "Outputs/Figures/Cluster_heatmap_2B"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# load data
metadata <- read_csv("Inputs/Metadata_discovery_WD_DD.csv")
exprs_data <- readRDS("Inputs/Batch_corrected_data_discovery_WD_DD.rds")
clusters <- read.csv("Outputs/Files/Consensus_clustering/Consensus_clustering.csv", row.names = 1)
degs <- read.csv("Outputs/Files/Limma/DEGs_6cl.csv")


metadata <- left_join(metadata,clusters, by = "sampleID")
colnames(metadata)[5] <- "cluster"

all_upregulated_genes <-  as.data.frame(degs)
all_upregulated_genes <- all_upregulated_genes %>%
  mutate(cluster = as.factor(cluster)) %>%
  mutate(cluster = case_when(
    cluster == "Cluster_1" ~ "A",
    cluster == "Cluster_2" ~ "B",
    cluster == "Cluster_3" ~ "C",
    cluster == "Cluster_4" ~ "D",
    cluster == "Cluster_5" ~ "E",
    cluster == "Cluster_6" ~ "F",
    TRUE ~ as.character(cluster)
  ))

#### complexpheatmap
metadata <- metadata[order(metadata$cluster), ]
all_upregulated_genes <- all_upregulated_genes[order(all_upregulated_genes$cluster),]
exprs_data_ordered <- exprs_data[, metadata$sampleID]
heatmap_data <- exprs_data_ordered[all_upregulated_genes$gene, ]

#  column annotation
annotation_col <- data.frame(
  Dataset    = factor(metadata$dataset),
  Differentiation = factor(metadata$subtype),
  Molecular_cluster = factor(metadata$cluster) )
rownames(annotation_col) <- colnames(heatmap_data)

# color lists
color_list <- c("#80B1D2","#FAF5B5","#BCB9D8","#F18072","#B3D46B","#F7CBDF","#8CD0C3","#F9B063","#D7D7D5","#BA7FB5" )
dataset_colors <- setNames(color_list[1:length(unique(metadata$dataset))],unique(metadata$dataset))

cluster_colors1 <- c(
  "1" = "#A5C1A4",
  "2" = "#827BB3",
  "3" = "#DF9677",
  "4" = "#1F70B5",
  "5" = "#F7CBDF",
  "6" = "#FDEDB9"
)

differentiation_colors <- c("DDLPS" = "#F18072", "WDLPS" = "#80B1D2", "WD/DDLPS" = "#F7CBDF")

ha_col <- HeatmapAnnotation(
  Dataset    = annotation_col$Dataset,
  Differentiation = annotation_col$Differentiation,
  Molecular_cluster = annotation_col$Molecular_cluster,
  annotation_name_side = "right",
  annotation_name_gp = gpar(fontsize = 10, fontface="bold"),
  col = list(
    Dataset    = dataset_colors,
    Differentiation = differentiation_colors,
    Molecular_cluster = cluster_colors1
  )
)

mat_zscore <- t(scale(t(heatmap_data)))

col_fun <- colorRamp2(c(-2, 0, 2), c("#80B1D2", "white", "darkred"))

out_file <- file.path(output_dir, "Cluster_heatmap.pdf")
# Open PDF device
pdf(out_file, width = 8, height = 6)

ht <- Heatmap(
  mat_zscore,
  name = "log2(CPM)\nZ-score",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  top_annotation = ha_col,
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 10, fontface = "bold"),
    legend_direction = "vertical"
  )
)

draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")

dev.off()

