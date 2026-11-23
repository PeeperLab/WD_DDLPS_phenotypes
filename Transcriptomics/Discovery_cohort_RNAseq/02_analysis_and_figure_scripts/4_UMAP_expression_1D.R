# This script creates a UMAP of the transcriptomic discovery cohort colored by the expression of a single gene.

library(ggplot2)
library(readr)
library(umap)

# set working directory
set.seed(1507)
## ---- Working directory ----------------------------------------------------
## Before running, set the working directory to the "Transcriptomics/Discovery_cohort_RNAseq" folder of
## this repository (the folder that holds Inputs/ and Outputs/), e.g.
## setwd("/path/to/WD_DDLPS_phenotypes/Transcriptomics/Discovery_cohort_RNAseq")

output_dir <- "Outputs/Figures/UMAP_expression_1D"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# load data
metadata <- read_csv("Inputs/Metadata_discovery_all.csv")
data <- readRDS("Inputs/Batch_corrected_data_discovery_all.rds")

gene_vars <- apply(data, 1, var)
top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:5000]
data_hvg <- data[top_genes, ]

### umap after batch effect correction
set.seed(1507)
# set UMAP parameters
umap_config <- umap.defaults
umap_config$n_neighbors <- 15  # adjust this parameter if necessary
umap_config$min_dist <- 0.8  # adjust this parameter if necessary

# UMAP
umap_result <- umap(t(data_hvg), config = umap_config)

umap_data <- as.data.frame(umap_result$layout)
colnames(umap_data) <- c("UMAP1", "UMAP2")

# match metadata rows to the samples (columns) of the expression matrix by sample ID
metadata <- metadata[match(colnames(data_hvg), metadata$sampleID), ]
stopifnot(identical(metadata$sampleID, colnames(data_hvg)))
umap_data$subtype <- metadata$subtype
umap_data$dataset <- metadata$dataset
umap_data$sampleID <- metadata$sampleID
umap_data$dataset <- as.factor(umap_data$dataset)

#plot single gene expression on UMAP
scaled_data <- scale(t(data))
scaled_data <- t(scaled_data)

plot_umap_by_gene <- function(gene, scaled_data, umap_data,
                              panel_width_in = 5.5, panel_height_in = 4) {
  
  # Check for gene existence
  if (!gene %in% rownames(scaled_data)) {
    stop(paste("Gene", gene, "not found in combat_data."))
  }
  
  # Ensure sample alignment
  if (!all(colnames(scaled_data) == rownames(umap_data))) {
    stop("Sample names in combat_data and umap_data do not match. Check rownames of umap_data.")
  }
  
  # Create expression column
  umap_data$expression <- as.numeric(scaled_data[gene, ])
  
  # Create plot
  p <- ggplot(umap_data, aes(x = UMAP1, y = -UMAP2, color = expression)) +
    geom_point(size = 3.8, alpha = 0.8, shape = 16) +
    coord_fixed() +
    scale_color_gradientn(colors = c("#80B1D2", "white", "darkred"),
                          name = paste0(gene, "\nexpression")) +
    theme_minimal() +
    labs(title = paste("UMAP Colored by", gene, "Expression")) +
    theme(
      text = element_text(size = 16),
      aspect.ratio = 1,
      axis.title = element_text(size = 16),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white"),
      panel.grid = element_blank())
  
  # Save to file
  ggsave(
    filename = file.path(output_dir, paste0("UMAP_", gene, "_plot.pdf")),
    plot = p,
    width = panel_width_in,
    height = panel_height_in,
    units = "in"
  )
  
  p
}

# Then call the function
# Fig. 1D
for (gene in c("MDM2", "PPARG", "TOP2A")) {
  plot_umap_by_gene(gene, scaled_data, umap_data)
}
