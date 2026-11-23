# This script re-plots the ORA results in bar graphs. 

library(tidyverse)
library(stringr)
library(ggplot2)

# set working directory
set.seed(1507)
## ---- Working directory ----------------------------------------------------
## Before running, set the working directory to the "Transcriptomics/Discovery_cohort_RNAseq" folder of
## this repository (the folder that holds Inputs/ and Outputs/), e.g.
## setwd("/path/to/WD_DDLPS_phenotypes/Transcriptomics/Discovery_cohort_RNAseq")

output_dir <- "Outputs/Figures/ORA/Re-plotted"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Find ORA CSV files
ora_csv_dir <- "Outputs/Files/ORA"
ora_files <- list.files(
  ora_csv_dir,
  pattern = "^rna_Cluster[0-9]+_enriched_pathways\\.csv$",
  full.names = TRUE
)

# Loop over clusters
for (file in ora_files) {
  
  # Extract cluster name (e.g. "Cluster3")
  cluster <- str_extract(basename(file), "Cluster[0-9]+")
  
  message("Plotting ORA for ", cluster)
  
  # Read ORA results
  ora_df <- read.csv(file, stringsAsFactors = FALSE)
  
  if (nrow(ora_df) == 0) {
    message("  No pathways found for ", cluster)
    next
  }
  

  # Prepare data
  ora_df <- ora_df %>%
    mutate(
      GeneRatio_num = as.numeric(sub("/.*", "", GeneRatio)) /
        as.numeric(sub(".*/", "", GeneRatio)),
      Description = str_replace_all(Description, "_", " ")
    ) %>%
    arrange(p.adjust)
  
  # Keep top pathways (optional)
  ora_df <- ora_df %>% slice_head(n = 15)
  
  ora_df <- ora_df %>%
    mutate(
      neglog10_padj = -log10(p.adjust)
    )
  
  midpoint_val <- -log10(0.05)
  
  upper_anchor <- quantile(
    ora_df$neglog10_padj,
    probs = 0.95,
    na.rm = TRUE
  )
  
  upper_anchor <- max(upper_anchor, midpoint_val + 0.5)
  
  # Plot
  p <- ggplot(
    ora_df %>% arrange(p.adjust),
    aes(
      x = GeneRatio_num,
      y = reorder(Description, -p.adjust),
      fill = neglog10_padj
    )
  ) +
    geom_bar(stat = "identity", color = "black", width = 0.8) +
    scale_fill_gradientn(
      colours = c("#80B1D2", "white", "darkred"),
      values = scales::rescale(c(0, midpoint_val, upper_anchor)),
      limits = c(0, upper_anchor),
      oob = scales::squish,
      name = expression(-log[10](adj.~p))
    ) +
    labs(
      title = paste(cluster, "– ORA enriched pathways"),
      x = "Gene Ratio",
      y = "Pathway"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.y = element_text(size = 9),
      plot.title = element_text(face = "bold", hjust = 0.5),
      panel.grid = element_blank(),                    
      panel.border = element_rect(color = "black",      
                                  fill = NA, 
                                  size = 1),
      panel.background = element_blank()             
    )
  
  # Save plot
  out_file <- file.path(
    output_dir,
    paste0("rna_", cluster, "_ORA_barplot.pdf")
  )
  
  ggsave(
    filename = out_file,
    plot = p,
    width = 10,
    height = max(4, 0.35 * nrow(ora_df)),
    dpi = 300
  )
}