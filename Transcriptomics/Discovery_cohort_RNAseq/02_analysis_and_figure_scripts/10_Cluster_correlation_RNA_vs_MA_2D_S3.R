# This script calculates the correlation of relative expression levels between clusters in the transcriptomic discovery cohort and the microarray validation cohort.
# For the validation cohort, limma results have been loaded in the Inputs folder, these can also be obtained by running the DEA script in the validation cohort folder.

library(tidyverse)
library(ggrepel)
library(pheatmap)
library(tibble)

# set working directory
set.seed(1507)
## ---- Working directory ----------------------------------------------------
## Before running, set the working directory to the "Transcriptomics/Discovery_cohort_RNAseq" folder of
## this repository (the folder that holds Inputs/ and Outputs/), e.g.
## setwd("/path/to/WD_DDLPS_phenotypes/Transcriptomics/Discovery_cohort_RNAseq")

# Input / output paths
RNA_dir <- "Outputs/Files/Limma"
MA_dir  <- "./Inputs/DEA_validation_files"

output_dir_files   <- "Outputs/Files/Cluster_correlations/RNA_vs_MA_cluster_correlations"
output_dir_figures <- "Outputs/Figures/Cluster_correlations_2D_S3"

dir.create(output_dir_files, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir_figures, recursive = TRUE, showWarnings = FALSE)

# Load DEA results
RNA_files <- list.files(RNA_dir, pattern = "Cluster_\\d+_DEGs", full.names = TRUE)
MA_files  <- list.files(MA_dir,  pattern = "Cluster_\\d+_DEGs", full.names = TRUE)

RNA_list <- lapply(RNA_files, read.csv, row.names = 1)
MA_list  <- lapply(MA_files,  read.csv, row.names = 1)

names(RNA_list) <- paste0("Cluster ", seq_along(RNA_list))
names(MA_list)  <- paste0("Cluster ", seq_along(MA_list))

# Compute correlations
correlation_results <- list()

for (i in seq_along(RNA_list)) {
  for (j in seq_along(MA_list)) {
    
    RNA <- RNA_list[[i]] %>%
      mutate(gene = toupper(rownames(.)))
    
    MA <- MA_list[[j]] %>%
      mutate(gene = toupper(rownames(.)))
    
    merged <- inner_join(RNA, MA, by = "gene",
                         suffix = c("_RNA", "_MA"))
    
    # Skip weak comparisons
    if (nrow(merged) < 20) next
    
    # Correlation
    cor_test <- cor.test(
      merged$logFC_RNA,
      merged$logFC_MA,
      method = "pearson"
    )
    
    r_value <- unname(cor_test$estimate)
    p_value <- cor_test$p.value
    
    # Top genes for labeling
    top_rna <- merged %>% slice_max(abs(t_RNA), n = 10)
    top_ma  <- merged %>% slice_max(abs(t_MA),  n = 10)
    
    top_genes <- bind_rows(top_rna, top_ma) %>%
      distinct(gene, .keep_all = TRUE)
    
    # Scatter plot
    p <- ggplot(merged, aes(x = logFC_RNA, y = logFC_MA)) +
      geom_point(alpha = 0.7, color = "#80B1D2", shape = 16) +
      geom_smooth(method = "lm", color = "darkred", se = FALSE) +
      geom_hline(yintercept = 0, color = "black", linewidth = 0.2) +
      geom_vline(xintercept = 0, color = "black", linewidth = 0.2) +
      theme_minimal() +
      labs(
        title = paste0(
          names(RNA_list)[i], " (Disc.) vs ",
          names(MA_list)[j], " (Val.)"
        ),
        x = "Log2(fold change) (Discovery)",
        y = "Log2(fold change) (Validation)"
      ) +
      annotate(
        "text",
        x = min(merged$logFC_RNA),
        y = max(merged$logFC_MA),
        label = paste0(
          "R = ", round(r_value, 3),
          "\nP = ", format.pval(p_value, digits = 3, eps = 0.001)
        ),
        hjust = 0,
        vjust = 1,
        size = 6,
        color = "black"
      ) +
      theme(
        panel.border     = element_rect(color = "black", fill = NA, size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title       = element_text(size = 17, face = "bold"),
        axis.title       = element_text(size = 15),
        axis.text        = element_text(size = 12)
      )
    
    # Save scatter plot
    plot_file <- file.path(
      output_dir_figures,
      paste0(
        names(RNA_list)[i], "_vs_",
        names(MA_list)[j], "_scatter.pdf"
      )
    )
    
    ggsave(plot_file, p, width = 4.5, height = 4.5)
    
    # Save stats
    correlation_results[[length(correlation_results) + 1]] <-
      data.frame(
        RNA_cluster = names(RNA_list)[i],
        MA_cluster  = names(MA_list)[j],
        R_value     = r_value,
        P_value     = p_value,
        N_genes     = nrow(merged)
      )
 }
}

correlation_results <- bind_rows(correlation_results)


write.csv(
  correlation_results,
  file.path(output_dir_files, "RNA_MA_cluster_correlations.csv"),
  row.names = FALSE
)

# Build correlation matrix
cor_matrix <- correlation_results %>%
  dplyr::select(RNA_cluster, MA_cluster, R_value) %>%
  tidyr::pivot_wider(names_from = MA_cluster, values_from = R_value) %>%
  tibble::column_to_rownames("RNA_cluster") %>%
  as.matrix()

RNA_order <- c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5", "Cluster 6")
MA_order  <- c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5", "Cluster 6")

cor_matrix_ordered <- cor_matrix[RNA_order, MA_order]

# Heatmap 1: ggplot (fixed order)
cor_long <- cor_matrix_ordered %>%
  as.data.frame() %>%
  rownames_to_column("RNA") %>%
  pivot_longer(-RNA, names_to = "MA", values_to = "R_value") %>%
  mutate(
    RNA = factor(RNA, levels = RNA_order),
    MA  = factor(MA,  levels = MA_order)
  )

heatmap_gg <- ggplot(cor_long, aes(RNA, MA, fill = R_value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "#80B1D2",
    mid = "white",
    high = "darkred",
    midpoint = 0,
    limits = c(-1, 1)
  ) +
  coord_fixed() +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    panel.grid = element_blank()
  ) +
  labs(
    title = "RNA–Microarray cluster correlation",
    x = "RNA clusters",
    y = "Microarray clusters",
    fill = "Pearson R"
  )

heatmap_gg

ggsave(
  file.path(output_dir_figures, "Correlation_heatmap_RNA_MA_ggplot.pdf"),
  heatmap_gg,
  width = 7, height = 6
)

