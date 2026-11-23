# This script creates stacks barplots that visualize the distribution of LPS subtypes accross individual datasets within the
# transcriptomic validation cohort. 

library(readr)
library(dplyr)
library(ggplot2)

# set working directory
set.seed(1507)
## ---- Working directory ----------------------------------------------------
## Before running, set the working directory to the "Transcriptomics/Validation_cohort_microarray" folder of
## this repository (the folder that holds Inputs/ and Outputs/), e.g.
## setwd("/path/to/WD_DDLPS_phenotypes/Transcriptomics/Validation_cohort_microarray")

output_dir <- "Outputs/Figures/Barplots_datasets_fig_S1B"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

## expression data
metadata <- read_csv("Inputs/Metadata_validation_all.csv")

# Proportional stacked barplot
barplot_data <- metadata %>%
  group_by(dataset, subtype) %>%
  summarise(N = n(), .groups = "drop")

barplot_data <- barplot_data %>%
  mutate(
    subtype = factor(
      subtype,
      levels = c("Fat_tissue", "WDLPS","DDLPS", "WD/DDLPS", "MLPS", "PLPS" )
    )
  )

stacked_barplot <- ggplot(
  barplot_data,
  aes(
    x = dataset,
    y = N,
    fill = subtype
  )
) +
  geom_col(position = "fill", size = 0.2, color = "black") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = c(
    "Fat_tissue" = "#FAF5B5",
    "WDLPS" = "#80B1D2",
    "DDLPS" = "#F18072",
    "WD/DDLPS" = "#F7CBDF",
    "MLPS" = "#4F7942",
    "PLPS" = "#060771"
  ), breaks = c("Fat_tissue", "WDLPS", "DDLPS", "WD/DDLPS", "MLPS", "PLPS")) +
  labs(
    x = "Dataset",
    y = "Fraction of samples",
    fill = "Subtype"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

stacked_barplot

ggsave(
  filename = file.path(output_dir, "Barplot_datasets_validation.pdf"),
  plot = stacked_barplot,
  width = 6,
  height = 5,
  units = "in"
)