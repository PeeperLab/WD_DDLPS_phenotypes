# This script performes correlation analysis between different clusters in the transcriptomic, and proteomic datasets. 
# The transcriptomic data is pre-loaded in the inputs folder of the proteomics, but can also be created by running the script in the transcriptomic dataset folders.

library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)

set.seed(1507)
## ---- Working directory ----------------------------------------------------
## Before running, set the working directory to the "Proteomics" folder of
## this repository (the folder that holds Inputs/ and Outputs/), e.g.
## setwd("/path/to/WD_DDLPS_phenotypes/Proteomics")

output_dir <- "Outputs/Files/Cluster_correlation"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

output_dir_2 <- "Outputs/Figures/Cluster_correlation_fig3B"
dir.create(output_dir_2, recursive = TRUE, showWarnings = FALSE)

# Load RNA datasets:
RNA_1 <- read.csv("Inputs/Transcriptomic_files/Cluster_1_DEGs_all.csv", stringsAsFactors = FALSE)
colnames(RNA_1)[colnames(RNA_1) == "X"] <- "gene"
RNA_2 <- read.csv("Inputs/Transcriptomic_files/Cluster_2_DEGs_all.csv", stringsAsFactors = FALSE)
colnames(RNA_2)[colnames(RNA_2) == "X"] <- "gene"
RNA_3 <- read.csv("Inputs/Transcriptomic_files/Cluster_3_DEGs_all.csv", stringsAsFactors = FALSE)
colnames(RNA_3)[colnames(RNA_3) == "X"] <- "gene"
RNA_4 <- read.csv("Inputs/Transcriptomic_files/Cluster_4_DEGs_all.csv", stringsAsFactors = FALSE)
colnames(RNA_4)[colnames(RNA_4) == "X"] <- "gene"
RNA_5 <- read.csv("Inputs/Transcriptomic_files/Cluster_5_DEGs_all.csv", stringsAsFactors = FALSE)
colnames(RNA_5)[colnames(RNA_5) == "X"] <- "gene"
RNA_6 <- read.csv("Inputs/Transcriptomic_files/Cluster_6_DEGs_all.csv", stringsAsFactors = FALSE)
colnames(RNA_6)[colnames(RNA_6) == "X"] <- "gene"

#protein datasets
Prot_1 <- read.csv("Outputs/Files/Limma/Cluster_1_DE_proteins_6cl.csv", stringsAsFactors = FALSE)
colnames(Prot_1)[colnames(Prot_1) == "X"] <- "gene"
Prot_2 <- read.csv("Outputs/Files/Limma/Cluster_2_DE_proteins_6cl.csv", stringsAsFactors = FALSE)
colnames(Prot_2)[colnames(Prot_2) == "X"] <- "gene"
Prot_3 <- read.csv("Outputs/Files/Limma/Cluster_3_DE_proteins_6cl.csv", stringsAsFactors = FALSE)
colnames(Prot_3)[colnames(Prot_3) == "X"] <- "gene"
Prot_4 <- read.csv("Outputs/Files/Limma/Cluster_4_DE_proteins_6cl.csv", stringsAsFactors = FALSE)
colnames(Prot_4)[colnames(Prot_4) == "X"] <- "gene"
Prot_5 <- read.csv("Outputs/Files/Limma/Cluster_5_DE_proteins_6cl.csv", stringsAsFactors = FALSE)
colnames(Prot_5)[colnames(Prot_5) == "X"] <- "gene"
Prot_6 <- read.csv("Outputs/Files/Limma/Cluster_6_DE_proteins_6cl.csv", stringsAsFactors = FALSE)
colnames(Prot_6)[colnames(Prot_6) == "X"] <- "gene"


# # List of RNA and protein datasets (6 RNA clusters + 6 protein clusters)
RNA_list <- list(RNA_1, RNA_2, RNA_3, RNA_4, RNA_5, RNA_6)
protein_list <- list(Prot_1, Prot_2, Prot_3, Prot_4, Prot_5, Prot_6)

# # Extract names for labeling (6 RNA clusters + 6 protein clusters)
RNA_names <- c("Cluster 1 (RNA)", "Cluster 2 (RNA)", "Cluster 3 (RNA)", "Cluster 4 (RNA)", "Cluster 5 (RNA)", "Cluster 6 (RNA)")
protein_names <- c("Cluster A (Protein)", "Cluster B (Protein)", "Cluster C (Protein)", "Cluster D (Protein)", "Cluster E (Protein)", "Cluster F (Protein)")

# Dataframe to store correlation results
correlation_results <- data.frame(RNA = character(), Protein = character(), R_value = numeric(), P_value = numeric(), stringsAsFactors = FALSE)

# Loop through RNA and protein datasets
for (i in seq_along(RNA_list)) {
  for (j in seq_along(protein_list)) {
    
    # Extract RNA & protein dataframe
    RNA <- RNA_list[[i]]
    protein <- protein_list[[j]]
    
    # Merge dataframes based on "gene" column
    merged_df <- inner_join(RNA, protein, by = "gene", suffix = c("_RNA", "_Protein"))
    
    # Perform correlation test
    cor_test <- cor.test(merged_df$logFC_RNA, merged_df$logFC_Protein, method = "pearson")
    r_value <- cor_test$estimate
    p_value <- cor_test$p.value
    
    # Identify top 10 genes with highest t-values in RNA and Protein
    top_rna_genes <- merged_df %>% top_n(10, t_RNA)  # Top 10 highest RNA t-values
    top_protein_genes <- merged_df %>% top_n(10, t_Protein)  # Top 10 highest Protein t-values
    top_genes <- bind_rows(top_rna_genes, top_protein_genes) %>% distinct(gene, .keep_all = TRUE)
    
    # Generate scatter plot
    plot <- ggplot(merged_df, aes(x = t_RNA, y = t_Protein)) +
      geom_point(alpha = 0.7, color = "#80B1D2", shape = 16) +
      #geom_text_repel(data = top_genes, aes(label = gene), vjust = -1, hjust = 1, size = 3, color = "darkred") +
      geom_smooth(method = "lm", color = "darkred", se = FALSE) +
      geom_hline(yintercept = 0, color = "black", linewidth = 0.2) +
      geom_vline(xintercept = 0, color = "black", linewidth = 0.2) +
      theme_minimal() +
      labs(
        title = paste(RNA_names[i], "vs", protein_names[j]),
        x = paste("Log2(fold change) ", RNA_names[i], sep = ""),
        y = paste("Log2(fold change) ", protein_names[j], sep = "")
      ) +
      annotate("text", x = min(merged_df$t_RNA), y = max(merged_df$t_Protein), 
               label = paste0("R = ", round(r_value, 3), "\nP = ", format.pval(p_value, digits = 3, eps = 0.001)), 
               hjust = 0, vjust = 1, size = 5, color = "black") +
      theme(
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title       = element_text(size = 14, face = "bold"),
        axis.title       = element_text(size = 12),
        axis.text        = element_text(size = 11)
      )
    
    # Save plot
    plot_filename <- paste0(output_dir, "/", RNA_names[i], "_vs_", protein_names[j], ".pdf")
    ggsave(plot_filename, plot, width = 4, height = 4, dpi = 300)
    
    # Save correlation results
    correlation_results <- rbind(correlation_results, data.frame(RNA = RNA_names[i], Protein = protein_names[j], R_value = r_value, P_value = p_value))
  }
}

# Save correlation results as CSV
write.csv(correlation_results, file = paste0(output_dir, "/correlation_results.csv"), row.names = FALSE)

# Load correlation results
# correlation_results <- read.csv("./Outputs/Files/Cluster_correlation/correlation_results.csv", stringsAsFactors = FALSE)

# Aggregate duplicates by averaging R-values
correlation_results <- correlation_results %>%
  group_by(RNA, Protein) %>%
  summarise(R_value = mean(R_value, na.rm = TRUE), .groups = "drop")  # Handle duplicates

# Convert data into a proper matrix format
cor_matrix <- correlation_results %>%
  pivot_wider(names_from = Protein, values_from = R_value) %>%
  column_to_rownames(var = "RNA")  # Set RNA as row names


# Define custom order
custom_order_RNA <- c(
  "Cluster 1 (RNA)", "Cluster 2 (RNA)", "Cluster 3 (RNA)",
  "Cluster 4 (RNA)", "Cluster 5 (RNA)", "Cluster 6 (RNA)"
)

custom_order_Protein <- c(
  "Cluster A (Protein)", "Cluster B (Protein)", "Cluster C (Protein)",
  "Cluster D (Protein)", "Cluster E (Protein)", "Cluster F (Protein)"
)

# Convert to long format and apply custom order
cor_long <- cor_matrix %>%
  rownames_to_column(var = "RNA") %>%
  pivot_longer(-RNA, names_to = "Protein", values_to = "R_value") %>%
  mutate(
    RNA = factor(RNA, levels = custom_order_RNA),
    Protein = factor(Protein, levels = custom_order_Protein)
  )

# Create heatmap with ggplot2
heatmap1 <- ggplot(cor_long, aes(x = Protein, y = RNA, fill = R_value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#80B1D2", mid = "white", high = "darkred", midpoint = 0) +
  theme_minimal() +
  labs(title = "Correlation Heatmap (RNA vs Protein Clusters)",
       x = "Protein Clusters",
       y = "RNA Clusters",
       fill = "R-value") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

heatmap1

# Save plot
ggsave(
  filename = file.path(output_dir_2, "Correlation_heatmap_RNA_protein.pdf"),
  plot = heatmap1,
  width = 8,
  height = 6,
  dpi = 300
)
