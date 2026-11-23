# GSEA of the Hallmark interferon-gamma response and a TLS signature in discovery cluster 1
# versus all other clusters (Fig. S5D).

library(tidyverse)
library(clusterProfiler)
library(msigdbr)
library(fgsea)

# set working directory
set.seed(1507)
## ---- Working directory ----------------------------------------------------
## Before running, set the working directory to the "Transcriptomics/Discovery_cohort_RNAseq" folder of
## this repository (the folder that holds Inputs/ and Outputs/), e.g.
## setwd("/path/to/WD_DDLPS_phenotypes/Transcriptomics/Discovery_cohort_RNAseq")

output_dir <- "Outputs/Figures/GSEA_cluster_1_S5D"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

output_dir_2 <- "Outputs/Files/GSEA_cluster_1_S5D"
if (!dir.exists(output_dir_2)) dir.create(output_dir_2, recursive = TRUE)

# Load cluster DEA results
cluster_files <- list.files(
  "Outputs/Files/Limma", 
  pattern = "Cluster_\\d+_DEGs_all.csv", full.names = TRUE
)

dea_results_list <- lapply(cluster_files, read.csv, row.names = 1)
names(dea_results_list) <- paste0("Cluster", seq_along(dea_results_list))

# Gene sets: MSigDB Hallmark collection and a custom tertiary lymphoid structure (TLS)
# signature (Sautes-Fridman et al., 2019). Only the gene sets in target_genesets are tested.
msigdb_gene_sets <- msigdbr(species = "Homo sapiens", category = "H") %>%
  select(gs_name, gene_symbol) %>%
  mutate(gene_symbol = toupper(gene_symbol))

custom_signatures <- read.csv("Inputs/GSEA_custom_signatures/TLS_Sautes-Fridman_2019.csv") %>%
  mutate(gene_symbol = toupper(gene_symbol))

gene_sets <- bind_rows(msigdb_gene_sets, custom_signatures)

# Select gene sets to test
target_genesets <- c(
  # MSigDB gene sets
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  # Custom signatures
  "TLS_Sautes-Fridman_2019"
)

cat("Testing gene sets:", paste(target_genesets, collapse = ", "), "\n\n")

# RUN GSEA
cluster_name <- "Cluster1"
results <- dea_results_list[[cluster_name]]

# Filter for significant genes
results_sig <- results %>%
  filter(!is.na(adj.P.Val), adj.P.Val < 0.01)

if(nrow(results_sig) < 10){
  cat("Too few significant genes for", cluster_name, "\n\n")
  stop("Insufficient significant genes for GSEA")
}

# Prepare gene list
gene_list <- results_sig$logFC
names(gene_list) <- toupper(rownames(results_sig))
gene_list <- gene_list[!is.na(gene_list)]
gene_list <- gene_list[!duplicated(names(gene_list))]
gene_list <- sort(gene_list, decreasing = TRUE)

cat("Gene list size:", length(gene_list), "\n\n")

# Filter gene sets
gene_sets_filtered <- gene_sets %>%
  filter(gs_name %in% target_genesets & gene_symbol %in% names(gene_list))

valid_genesets <- gene_sets_filtered %>%
  group_by(gs_name) %>%
  summarise(n_genes = n(), .groups = "drop") %>%
  filter(n_genes >= 1) %>%
  pull(gs_name)

gene_sets_filtered <- gene_sets_filtered %>% 
  filter(gs_name %in% valid_genesets)


if(length(valid_genesets) == 0){
  stop("No valid gene sets with sufficient genes")
}

# Run GSEA
gsea_res <- GSEA(
  gene_list,
  exponent = 1,
  minGSSize = 1,
  maxGSSize = 500,
  pvalueCutoff = 1,
  TERM2GENE = gene_sets_filtered,
  TERM2NAME = NA,
  verbose = TRUE,
  seed = FALSE
)

res_df <- as.data.frame(gsea_res)

if(nrow(res_df) == 0){
  cat("No results from GSEA\n")
  stop("No enrichment found")
}

# Function to extract curve and ticks
getEnrichmentData <- function(pathway_genes, stats){
  pd <- plotEnrichmentData(pathway = pathway_genes, stats = stats, gseaParam = 1)
  return(pd)
}

# Prepare data for combined plot
plot_list <- lapply(valid_genesets, function(gs){
  genes <- gene_sets_filtered$gene_symbol[gene_sets_filtered$gs_name == gs]
  pd <- getEnrichmentData(genes, gene_list)
  pd$GeneSet <- gs
  return(pd)
})

# Combine curves and ticks
tick_spacing <- 0.05
floor_val <- 0

curve_df <- map_dfr(plot_list, ~ { .x$curve %>% mutate(GeneSet = .x$GeneSet) })
ticks_df <- map_dfr(seq_along(plot_list), function(i){
  pd <- plot_list[[i]]
  tick_y_top <- floor_val - (i-1)*tick_spacing
  tick_height <- tick_spacing * 0.8
  pd$ticks %>% mutate(
    y_start = tick_y_top,
    y_end = tick_y_top - tick_height,
    GeneSet = pd$GeneSet
  )
})

# Prepare significance labels
sig_labels <- res_df %>%
  filter(ID %in% valid_genesets) %>%
  mutate(asterisk = case_when(
    p.adjust < 0.001 ~ "***",
    p.adjust < 0.01  ~ "**",
    p.adjust < 0.05  ~ "*",
    TRUE ~ ""
  )) %>%
  mutate(GeneSetLabel = paste(ID, asterisk, sep=" ")) %>%
  select(ID, GeneSetLabel)

curve_df <- curve_df %>%
  left_join(sig_labels, by=c("GeneSet"="ID"))
ticks_df <- ticks_df %>%
  left_join(sig_labels, by=c("GeneSet"="ID"))

# Build plot
colors <- RColorBrewer::brewer.pal(max(3, length(valid_genesets)), "Set1")[1:length(valid_genesets)]

p <- ggplot() +
  geom_line(data=curve_df, aes(x=rank, y=ES, color=GeneSetLabel), size=1.2) +
  geom_segment(data=ticks_df, aes(x=rank, xend=rank, y=y_start, yend=y_end, color=GeneSetLabel), size=0.5) +
  geom_hline(yintercept=0, color="black") +
  labs(
    title=paste(cluster_name, "- Combined GSEA"),
    x="Rank",
    y="Enrichment Score (ES)",
    color="Gene Set"
  ) +
  theme_minimal() +
  theme(
    axis.text=element_text(size=10),
    axis.title=element_text(size=12, face="bold"),
    legend.text=element_text(size=10),
    legend.title=element_text(size=12, face="bold")
  ) +
  scale_color_manual(values=colors)

print(p)

ggsave(
  filename = file.path(output_dir, "GSEA_cluster_1.pdf"),
  plot = p,
  width = 12,
  height = 6,
  units = "in")

write.csv(res_df, file = file.path(output_dir_2, "GSEA_cluster_1.csv"), row.names = TRUE)
