# This script performs GSEA on a selection of genesets from msigdb on cluster 4.

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(clusterProfiler)
library(msigdbr)
library(fgsea)
library(RColorBrewer)
library(purrr)  # map_dfr()

# set working directory
set.seed(1507)
## ---- Working directory ----------------------------------------------------
## Before running, set the working directory to the "Transcriptomics/Discovery_cohort_RNAseq" folder of
## this repository (the folder that holds Inputs/ and Outputs/), e.g.
## setwd("/path/to/WD_DDLPS_phenotypes/Transcriptomics/Discovery_cohort_RNAseq")

# Load cluster DEA results
cluster_files <- list.files(
  "./Outputs/Files/Limma",
  pattern = "Cluster_\\d+_DEGs_all.csv",
  full.names = TRUE
)

# Create a named list of DEA results
dea_results_list <- lapply(cluster_files, read.csv, row.names = 1)
names(dea_results_list) <- paste0("Cluster", seq_along(dea_results_list))

# Gene sets: MSigDB C2 (Reactome, WikiPathways) and C5 GO:BP.
# Only the gene sets in target_genesets are tested.
gene_sets <- bind_rows(
  msigdbr(species = "Homo sapiens", category = "C2"),
  msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")
) %>%
  select(gs_name, gene_symbol) %>%
  mutate(gene_symbol = toupper(gene_symbol))

# Gene sets of interest
target_genesets <- c(
  "REACTOME_SIGNALING_BY_WNT",
  "GOBP_WNT_SIGNALING_PATHWAY",
  "WP_WNT_SIGNALING",
  "REACTOME_TCF_DEPENDENT_SIGNALING_IN_RESPONSE_TO_WNT",
  "GOBP_CANONICAL_WNT_SIGNALING_PATHWAY",
  "GOBP_STEM_CELL_DEVELOPMENT",
  "GOBP_EMBRYONIC_ORGAN_MORPHOGENESIS"
)

# Select cluster
cluster_name <- "Cluster4"
results <- dea_results_list[[cluster_name]]

# Filter for significant genes
results_sig <- results %>%
  filter(!is.na(adj.P.Val), adj.P.Val < 0.01)

if(nrow(results_sig) < 10){
  cat("Too few significant genes for", cluster_name, "\n\n")
  stop("Insufficient significant genes for GSEA")
}

# Prepare gene list (from significant genes only)
gene_list <- results_sig$logFC
names(gene_list) <- toupper(rownames(results_sig))
gene_list <- gene_list[!is.na(gene_list)]
gene_list <- gene_list[!duplicated(names(gene_list))]
gene_list <- sort(gene_list, decreasing = TRUE)

# Filter gene sets
gene_sets_filtered <- gene_sets %>%
  filter(gs_name %in% target_genesets & gene_symbol %in% names(gene_list))

valid_genesets <- gene_sets_filtered %>%
  group_by(gs_name) %>%
  summarise(n_genes = n(), .groups = "drop") %>%
  filter(n_genes >= 1) %>%
  pull(gs_name)

gene_sets_filtered <- gene_sets_filtered %>% filter(gs_name %in% valid_genesets)

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

# Combine curves and ticks into single dataframes
tick_spacing <- 0.05
floor_val <- 0  # baseline for first row

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

# Prepare significance labels for legend (asterisks)
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

# Update curve_df and ticks_df with labeled GeneSet
curve_df <- curve_df %>%
  left_join(sig_labels, by=c("GeneSet"="ID"))
ticks_df <- ticks_df %>%
  left_join(sig_labels, by=c("GeneSet"="ID"))

# Build combined ggplot
colors <- RColorBrewer::brewer.pal(length(valid_genesets), "Set1")

p <- ggplot() +
  geom_line(data=curve_df, aes(x=rank, y=ES, color=GeneSetLabel), size=1.2) +
  geom_segment(data=ticks_df, aes(x=rank, xend=rank, y=y_start, yend=y_end, color=GeneSetLabel), size=0.5) +
  geom_hline(yintercept=0, color="black") +
  labs(title=paste(cluster_name, "- Combined GSEA"),
       x="Rank",
       y="Enrichment Score (ES)",
       color="Gene Set") +
  theme_minimal() +
  theme(
    axis.text=element_text(size=10),
    axis.title=element_text(size=12, face="bold"),
    legend.text=element_text(size=10),
    legend.title=element_text(size=12, face="bold")
  ) +
  scale_color_manual(values=colors)

print(p)

# Save plot
gsea_figdir <- "Outputs/Figures/GSEA_cluster_4_S6A"
dir.create(gsea_figdir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(gsea_figdir, paste0(cluster_name, "_GSEA_cluster_4.pdf")),
       plot=p, width=12, height=6)

# Save GSEA results table
gsea_filedir <- "Outputs/Files/GSEA_cluster_4_S6A"
dir.create(gsea_filedir, recursive = TRUE, showWarnings = FALSE)
results_table <- res_df %>%
  select(ID, Description, setSize, enrichmentScore, NES, pvalue, p.adjust, qvalue, rank) %>%
  arrange(p.adjust)

write.csv(results_table,
          file.path(gsea_filedir, paste0(cluster_name, "_GSEA_cluster_4_results.csv")),
          row.names = FALSE)