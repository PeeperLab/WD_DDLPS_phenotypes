# This script creates a heatmap of the ssGSEA values per cluster. 

library(readr)
library(ComplexHeatmap)
library(circlize)
library(grid)

set.seed(1507)
## ---- Working directory ----------------------------------------------------
## Before running, set the working directory to the "Proteomics" folder of
## this repository (the folder that holds Inputs/ and Outputs/), e.g.
## setwd("/path/to/WD_DDLPS_phenotypes/Proteomics")

output_dir <- "Outputs/Figures/ssGSEA_fig_3C"
if (!dir.exists(output_dir)) {dir.create(output_dir, recursive = TRUE)}

## ssGSEA result pathway
file_path <- "./Outputs/Files/ssGSEA/gsva_results_RNA_sigs.csv"
metadata <- read.csv("Outputs/Files/Consensus_clustering/Consensus_clustering.csv", row.names = 1)

## ssGSEA data
ssGSEA_res <- read_csv(file_path)
ssGSEA_res <- as.data.frame(ssGSEA_res)
rownames(ssGSEA_res) <- ssGSEA_res$pathway
ssGSEA_res <- ssGSEA_res[,-1]
metadata$sampleID <- rownames(metadata)
ssGSEA_res <- ssGSEA_res[ ,metadata$sampleID]

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

desired_order <- c("1", "2", "3", "4", "5", "6")
sub_metadata <- metadata
sub_metadata$consensus_clusters <- factor(sub_metadata$consensus_clusters, levels = desired_order)
sub_metadata <- sub_metadata[order(sub_metadata$consensus_clusters), ]
ssGSEA_res_ordered <- ssGSEA_res[, sub_metadata$sampleID]

#### average the ssGSEA score for each subtype
subtype <- unique(sub_metadata$consensus_clusters)

mean_matrix <- sapply(subtype, function(st) {
  idx <- which(sub_metadata$consensus_clusters == st)
  rowMeans(ssGSEA_res_ordered[, idx, drop = FALSE])
})
rownames(mean_matrix) <- rownames(ssGSEA_res_ordered)
colnames(mean_matrix) <- subtype
mat_z <- t(scale(t(mean_matrix)))
col_fun2 <- colorRamp2(c(min(mat_z), 0, max(mat_z)),
                       c("#80B1D2",  "#FFFFFF", "darkred"))


out_file <- file.path(output_dir, "ssGSEA_heatmap.pdf")
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)

# Prepare rownames
pretty_rownames <- gsub("_", " ", rownames(mat_z))

# Open PDF device
pdf(out_file, width = 8, height = 6)

# Create heatmap object
ht2 <- Heatmap(
  mat_z,
  name               = "z-score",
  col                = col_fun2,
  cluster_rows       = FALSE,
  cluster_columns    = FALSE,
  show_row_names     = TRUE,
  show_column_names  = TRUE,
  row_names_gp       = gpar(fontsize = 12),
  row_labels         = pretty_rownames,
  row_names_max_width = unit(16, "cm"),
  row_names_side     = "right",
  row_names_rot      = 0,
  column_names_gp    = gpar(fontsize = 4, fontface = "bold"),
  column_names_rot   = -45,
  heatmap_legend_param = list(
    title            = "Z-score of Mean-ssSGEA score",
    legend_direction = "horizontal",
    legend_width     = unit(1.5, "cm"),
    grid_height      = unit(1, "mm"),
    grid_width       = unit(1, "mm"),
    title_gp         = gpar(fontsize = 4, fontface = "bold"),
    labels_gp        = gpar(fontsize = 4)
  )
)

draw(ht2,
     heatmap_legend_side    = "top",
     annotation_legend_side = "top")

dev.off()
