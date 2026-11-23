# This script calculates the correlation of relative expression levels between clusters in the transcriptomic validation cohort

library(tidyverse)
library(ggrepel)
library(pheatmap)
library(tibble)

# set working directory
set.seed(1507)
## ---- Working directory ----------------------------------------------------
## Before running, set the working directory to the "Transcriptomics/Validation_cohort_microarray" folder of
## this repository (the folder that holds Inputs/ and Outputs/), e.g.
## setwd("/path/to/WD_DDLPS_phenotypes/Transcriptomics/Validation_cohort_microarray")

# Input / output paths
RNA_dir <- "Outputs/Files/Limma"

output_dir_files   <- "Outputs/Files/Cluster_correlations/MA_vsMA_cluster_correlations"
output_dir_figures <- "Outputs/Figures/Cluster_correlations_S2F"

dir.create(output_dir_files, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir_figures, recursive = TRUE, showWarnings = FALSE)

# Load DEA results
RNA_files <- list.files(RNA_dir, pattern = "Cluster_\\d+_DEGs", full.names = TRUE)

RNA_list <- lapply(RNA_files, read.csv, row.names = 1)
names(RNA_list) <- paste0("Cluster", seq_along(RNA_list))


# Compute correlations
correlation_results <- list()

for (i in seq_along(RNA_list)) {
  for (j in seq_along(RNA_list)) {
    
    RNA_A <- RNA_list[[i]] %>%
      mutate(gene = toupper(rownames(.)))
    
    RNA_B <- RNA_list[[j]] %>%
      mutate(gene = toupper(rownames(.)))
    
    merged <- inner_join(
      RNA_A, RNA_B,
      by = "gene",
      suffix = c("_A", "_B")
    )
    
    # Skip weak comparisons
    if (nrow(merged) < 20) next
    
    # Correlation
    cor_test <- cor.test(
      merged$logFC_A,
      merged$logFC_B,
      method = "pearson"
    )
    
    r_value <- unname(cor_test$estimate)
    p_value <- cor_test$p.value
    
    # Scatter plot
    p <- ggplot(merged, aes(x = logFC_A, y = logFC_B)) +
      geom_point(alpha = 0.7, color = "#80B1D2", shape = 16) +
      geom_smooth(method = "lm", color = "darkred", se = FALSE) +
      geom_hline(yintercept = 0, color = "black", linewidth = 0.2) +
      geom_vline(xintercept = 0, color = "black", linewidth = 0.2) +
      theme_minimal() +
      labs(
        title = paste0(
          names(RNA_list)[i], " vs ",
          names(RNA_list)[j]
        ),
        x = paste0(names(RNA_list)[i], " log2FC"),
        y = paste0(names(RNA_list)[j], " log2FC")
      ) +
      annotate(
        "text",
        x = min(merged$logFC_A),
        y = max(merged$logFC_B),
        label = paste0(
          "R = ", round(r_value, 3),
          "\nP = ", format.pval(p_value, digits = 3, eps = 0.001)
        ),
        hjust = 0,
        vjust = 1,
        size = 5
      ) +
      theme(
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        panel.grid = element_blank()
      )
    
    ggsave(
      file.path(
        output_dir_figures,
        paste0(
          names(RNA_list)[i], "_vs_",
          names(RNA_list)[j], "_scatter.pdf"
        )
      ),
      p,
      width = 4,
      height = 4
    )
    
    # Store stats
    correlation_results[[length(correlation_results) + 1]] <-
      data.frame(
        Cluster_A = names(RNA_list)[i],
        Cluster_B = names(RNA_list)[j],
        R_value   = r_value,
        P_value   = p_value,
        N_genes   = nrow(merged)
      )
  }
}

correlation_results <- bind_rows(correlation_results)

write.csv(
  correlation_results,
  file.path(output_dir_files, "MA_MA_cluster_correlations.csv"),
  row.names = FALSE
)

# Build correlation matrix
cor_matrix <- correlation_results %>%
  select(Cluster_A, Cluster_B, R_value) %>%
  pivot_wider(names_from = Cluster_B, values_from = R_value) %>%
  column_to_rownames("Cluster_A") %>%
  as.matrix()

# Manual cluster ordering
RNA_order <- paste0("Cluster", 1:6)

cor_matrix_ordered <- cor_matrix[RNA_order, RNA_order]

# MASK upper triangle
cor_matrix_tri <- cor_matrix_ordered
cor_matrix_tri[upper.tri(cor_matrix_tri)] <- NA

# Convert to long format
cor_long <- cor_matrix_tri %>%
  as.data.frame() %>%
  rownames_to_column("Cluster_A") %>%
  pivot_longer(-Cluster_A, names_to = "Cluster_B", values_to = "R_value") %>%
  mutate(
    Cluster_A = factor(Cluster_A, levels = RNA_order),
    Cluster_B = factor(Cluster_B, levels = RNA_order)
  )

# Triangular heatmap
heatmap_gg <- ggplot(cor_long, aes(Cluster_A, Cluster_B, fill = R_value)) +
  geom_tile(color = "white", na.rm = TRUE) +
  scale_fill_gradient2(
    low = "#80B1D2",
    mid = "white",
    high = "darkred",
    midpoint = 0,
    limits = c(-1, 1),
    na.value = "white"
  ) +
  coord_fixed() +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    panel.grid = element_blank()
  ) +
  labs(
    title = "Microarray–microarray cluster correlation",
    x = "Cluster A",
    y = "Cluster B",
    fill = "Pearson R"
  )

heatmap_gg

ggsave(
  file.path(output_dir_figures, "Correlation_heatmap_MA_MA.pdf"),
  heatmap_gg,
  width = 7, height = 6
)
