# This script creates UMAPS of the transcriptomic discovery cohort colored by subtype and dataset, 
# and a UMAP with WD/DDLPS samples highlighted in color and non-WD/DDLPS samples in grey.

library(ggplot2)
library(readr)
library(umap)

# set working directory
set.seed(1507)
## ---- Working directory ----------------------------------------------------
## Before running, set the working directory to the "Transcriptomics/Discovery_cohort_RNAseq" folder of
## this repository (the folder that holds Inputs/ and Outputs/), e.g.
## setwd("/path/to/WD_DDLPS_phenotypes/Transcriptomics/Discovery_cohort_RNAseq")

output_dir <- "Outputs/Figures/UMAP_all_samples_1C_S1C"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

output_dir_2 <- "Outputs/Figures/UMAP_subset_highlight_2A"
if (!dir.exists(output_dir_2)) dir.create(output_dir_2, recursive = TRUE)

# load data
metadata <- read_csv("Inputs/Metadata_discovery_all.csv")
data <- readRDS("Inputs/Batch_corrected_data_discovery_all.rds")

# add subset metadata for highlighting WD/DDLPS samples in UMAP
metadata_subset <- read_csv("Inputs/Metadata_discovery_WD_DD.csv")

gene_vars <- apply(data, 1, var)
top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:5000]
data_hvg <- data[top_genes, ]

# set UMAP parameters
umap_config <- umap.defaults
umap_config$n_neighbors <- 15 
umap_config$min_dist <- 0.8 

# UMAP
umap_result <- umap(t(data_hvg), config = umap_config)
umap_data <- as.data.frame(umap_result$layout)
colnames(umap_data) <- c("UMAP1", "UMAP2")

# add metadata information to umap_data
# match metadata rows to the samples (columns) of the expression matrix by sample ID
metadata <- metadata[match(colnames(data_hvg), metadata$sampleID), ]
stopifnot(identical(metadata$sampleID, colnames(data_hvg)))
umap_data$subtype <- metadata$subtype
umap_data$dataset <- metadata$dataset
umap_data$sampleID <- metadata$sampleID
umap_data$dataset <- as.factor(umap_data$dataset)
umap_data$sampleID <- rownames(umap_data)

#plot UMAP colored by subtype
subtype_colors <- c("DDLPS" = "#F18072", "WDLPS" = "#80B1D2", "Fat_tissue" = "#FAF5B5", "PLPS" = "#060771", "MLPS" = "#4F7942", "WD/DDLPS" = "#F7CBDF")

umap <- ggplot(umap_data, aes(x = UMAP1, y = -UMAP2, color = subtype, label = NA)) +
  geom_point(size = 5.5, alpha = 0.8, shape = 16) +
  coord_fixed() +
  theme_minimal() +
  labs(title = "", color = "Sample type") +
  theme(legend.position = "right",
        aspect.ratio = 1,
        panel.background = element_rect(fill = "white",colour = NA),
        plot.background = element_rect(fill = "white",colour = NA),
        text = element_text(size = 16),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16, face = "bold"),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        panel.grid = element_blank())+
  scale_color_manual(values = subtype_colors
  ) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

umap

ggsave(
  filename = file.path(output_dir, "UMAP_all_samples_by_subtype_discovery.pdf"),
  plot = umap,
  width = 7.5,
  height = 6,
  units = "in"
)

#plot UMAP colored by dataset
umap2 <- ggplot(umap_data, aes(x = UMAP1, y = -UMAP2, color = dataset, label = NA)) +
  geom_point(size = 5.5, alpha = 0.8, shape = 16) +
  coord_fixed() +
  theme_minimal() +
  labs(title = "", color = "Sample type") +
  theme(legend.position = "right",
        aspect.ratio = 1,
        panel.background = element_rect(fill = "white",colour = NA),
        plot.background = element_rect(fill = "white",colour = NA),
        text = element_text(size = 16),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16, face = "bold"),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        panel.grid = element_blank())+
  scale_color_brewer(palette = "Set3" 
  ) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

umap2

ggsave(
  filename = file.path(output_dir, "UMAP_all_samples_by_dataset_discovery.pdf"),
  plot = umap2,
  width = 7.5,
  height = 6,
  units = "in"
)

################## UMAP with non-WD/DDLPS samples colored in grey ####################
umap3 <- ggplot() +
  # ---- background: non-subset samples in grey ----
geom_point(
  data = subset(umap_data, !(sampleID %in% metadata_subset$sampleID)),
  aes(x = UMAP1, y = -UMAP2),
  color = "grey80",
  size = 4,
  alpha = 0.4,
  shape = 16
) +

  # ---- foreground: subset samples colored by subtype ----
geom_point(
  data = subset(umap_data, sampleID %in% metadata_subset$sampleID),
  aes(x = UMAP1, y = -UMAP2, color = subtype),
  size = 5.5,
  alpha = 0.8,
  shape = 16
) +
  coord_fixed() +
  scale_color_manual(values = c(
    "DDLPS" = "#F18072",
    "WDLPS" = "#80B1D2",
    "Fat_tissue" = "#FAF5B5",
    "PLPS" = "#060771",
    "MLPS" = "#4F7942",
    "WD/DDLPS" = "#F7CBDF"
  )) +
  theme_minimal() +
  labs(color = "Sample type") +
  theme(
    legend.position = "right",
    aspect.ratio = 1,
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

umap3

ggsave(
  filename = file.path(output_dir_2, "UMAP_all_samples_subset_highlight_discovery.pdf"),
  plot = umap3,
  width = 7.5,
  height = 6,
  units = "in"
)
