# This scripts creates UMAPS using all samples for the proteomic dataset

library(ggplot2)
library(dplyr)
library(openxlsx)
library(umap)

set.seed(1507)
## ---- Working directory -
## Before running, set the working directory to the "Proteomics" folder of
## this repository (the folder that holds Inputs/ and Outputs/), e.g.
## setwd("/path/to/WD_DDLPS_phenotypes/Proteomics")

proteomic_data <- read.csv("Inputs/Proteomic_data_batch_corrected.csv", row.names = 1)
proteomic_metadata <- read.xlsx("Inputs/Proteomics_metadata.xlsx")

output_dir <- "Outputs/Figures/UMAP_all_samples_Fig_1H_S1FG"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Put ID as rownames and order
rownames(proteomic_metadata) <- proteomic_metadata[[1]] 

## UMAP all samples (Fig. 1H) 
umap_data <- as.matrix(proteomic_data)

umap_data <- t(umap_data)
protein_variance <- apply(umap_data, 1, var)
selected_proteins <- names(sort(protein_variance, decreasing = TRUE))[1:500]
umap_data <- umap_data[selected_proteins, ]
umap_data <- t(umap_data)

set.seed(1507) 
n_neighbors <- 25
min_dist <- 0.8

umap_result <- umap(umap_data, n_neighbors = n_neighbors, min_dist = min_dist)

# Extract UMAP coordinates
umap_coords <- as.data.frame(umap_result$layout)
colnames(umap_coords) <- c("UMAP1", "UMAP2")
umap_coords$SampleID <- rownames(umap_data)

# Visualize the UMAP after batch correction
# Add metadata matched by sample ID
meta_matched <- proteomic_metadata[umap_coords$SampleID, ]
stopifnot(identical(meta_matched$ID, umap_coords$SampleID))
umap_coords$class <- meta_matched$class
umap_coords$batch <- meta_matched$batch
umap_coords$pat <- meta_matched$pat
umap_coords$ID <- meta_matched$ID


# Plot using ggplot2
umap_all <- ggplot(umap_coords, aes(x = -UMAP1, y = UMAP2, color = class)) +
  geom_point(alpha = 0.8, size = 6, shape = 16) +
  coord_fixed() +
  theme_minimal() +
   scale_color_manual(values = c(
    "DD" = "#F18072",
    "WD" = "#80B1D2",
    "Healthy_Fat" = "#FAF5B5",
    "Lipoma" = "#F7CBDF")) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    aspect.ratio = 1,
    panel.grid = element_blank())

umap_all

#Save as PDF
ggsave(
  filename = file.path(output_dir, "UMAP_Fig_1H.pdf"),
  plot = umap_all,
  width = 8,
  height = 6,
  dpi = 300
)

# # # Highlight specific samples on the UMAP

highlight <- umap_coords[umap_coords$pat %in% c(15, 25, 38), ]
other <- umap_coords[!umap_coords$pat %in% c(15, 25, 38), ]

# Plot
umap_highlight <- ggplot() +
  geom_point(data = other, aes(x = -UMAP1, y = UMAP2, color = class), alpha = 0.2, size = 6, shape = 16) +
  geom_point(data = highlight, aes(x = -UMAP1, y = UMAP2, color = pat), alpha = 1, size = 6, shape = 16) +
  coord_fixed() +
  theme_minimal() +
  scale_color_manual(values = c(
    "DD" = "#F18072",
    "WD" = "#80B1D2",
    "Healthy_Fat" = "#FAF5B5",
    "Lipoma" = "#F7CBDF",
    "15" = "#7AC36A",
    "25" = "#4DBBD5",
    "38" = "#9E4DB5")) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    aspect.ratio = 1,
    panel.grid = element_blank())
umap_highlight

#Save as PDF
ggsave(
  filename = file.path(output_dir, "UMAP_Fig_S1G.pdf"),
  plot = umap_highlight,
  width = 8,
  height = 6,
  dpi = 300
)


# ###### UMAP with protein expression 
# Scale each column (protein) across samples
proteomic_data_scaled <- as.data.frame(scale(proteomic_data))
# Fig. S1F
for (protein_of_interest in c("CDK4", "ADIPOQ", "TOP2B")) {

  # Ensure scaled proteomic data includes sample IDs as a column
  proteomic_data_scaled$SampleID <- rownames(proteomic_data_scaled)

  # Extract expression of the selected protein
  expression_df <- proteomic_data_scaled %>%
    select(SampleID, all_of(protein_of_interest)) %>%
    setNames(c("SampleID", "expression"))

  # Merge with UMAP coordinates by SampleID
  umap_expr <- umap_coords %>%
    left_join(expression_df, by = "SampleID")

  # Create the plot
  p <- ggplot(umap_expr, aes(x = -UMAP1, y = UMAP2, color = expression)) +
    geom_point(size = 5, alpha = 0.8, shape = 16) +
    coord_fixed() +
    scale_color_gradientn(
      colors = c("#80B1D2", "white", "darkred"),
      name = paste0(protein_of_interest, "\nExpression (scaled)")
    ) +
    theme_minimal() +
    labs(
      title = paste("UMAP Colored by", protein_of_interest, "Expression"),
      x = "UMAP 1",
      y = "UMAP 2"
    ) +
    theme(
      text = element_text(size = 16),
      aspect.ratio = 1,
      axis.title = element_text(size = 16),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 12),
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      panel.background = element_rect(fill = "white"),
      panel.grid = element_blank(),
      plot.background = element_rect(fill = "white")
    )
  p

  # Save as PDF
  ggsave(
    filename = file.path(output_dir, paste0("UMAP_", protein_of_interest, "_expression_Fig_S1F.pdf")),
    plot = p,
    width = 7,
    height = 5,
    dpi = 300
  )
}
