# This script creates a treemap visualization of individual datasets and subtypes within the transcriptomic validation cohort. 

library(readr)
library(dplyr)
library(ggplot2)
library(treemapify)

# set working directory
set.seed(1507)
## ---- Working directory ----------------------------------------------------
## Before running, set the working directory to the "Transcriptomics/Validation_cohort_microarray" folder of
## this repository (the folder that holds Inputs/ and Outputs/), e.g.
## setwd("/path/to/WD_DDLPS_phenotypes/Transcriptomics/Validation_cohort_microarray")

output_dir <- "Outputs/Figures/Treemap_fig_1B"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

## expression data
metadata <- read_csv("Inputs/Metadata_validation_all.csv")

treemap_data <- metadata %>%
  group_by(subtype, dataset) %>%
  summarise(N = n(), .groups = "drop")

# Aggregate metadata
treemap_data <- metadata %>%
  group_by(subtype, dataset) %>%
  summarise(N = n(), .groups = "drop")


# Treemap with hierarchical grouping (Subtype → Dataset)
treemap <- ggplot(
  treemap_data,
  aes(
    area = N,
    fill = subtype,
    subgroup = subtype,
    subgroup2 = dataset,
    label = paste0(dataset, "\nN = ", N)
  )
) +
  # big colored subtype blocks
  geom_treemap(color = "white", size = 2,) +
  
  # outline subtype groups
  geom_treemap_subgroup_border(color = "black", size = 2) +
  
  # outline dataset blocks within each subtype
  geom_treemap_subgroup2_border(color = "black", size = 2) +
  
  # dataset labels inside the subtype blocks
  geom_treemap_text(
    place = "centre",
    grow = TRUE,       # allows the text to scale to fill the block
    reflow = FALSE,    # prevents automatic line breaks
    colour = "black",
    min.size = 3,     
  ) +
  
  scale_fill_manual(values = c(
    "DDLPS" = "#E5988F",
    "WDLPS" = "#8CA8D9",
    "Fat_tissue" = "#FAF5B5",
    "PLPS" = "#8D6AAE",
    "MLPS" = "#9BCB7A",
    "WD/DDLPS" = "#F7CBDF"
  )) +
  coord_transform(y = "reverse") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    panel.grid = element_blank()
  )

treemap

ggsave(
  filename = file.path(output_dir, "Treemap_validation.pdf"),
  plot = treemap,
  width = 6,
  height = 8,
  units = "in"
)