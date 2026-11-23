# This script performs overrepresentation analysis for each cluster. It uses an external script to perform ORA, located in the code folder.
# Results are saved as tables; the figure panels (Fig. 4B, 5B) are drawn by 20_ORA_replotting_4B_5B.R.

library(tidyverse)
library(clusterProfiler)
library(msigdbr)
library(limma)

# set working directory
set.seed(1507)
## ---- Working directory ----------------------------------------------------
## Before running, set the working directory to the "Transcriptomics/Discovery_cohort_RNAseq" folder of
## this repository (the folder that holds Inputs/ and Outputs/), e.g.
## setwd("/path/to/WD_DDLPS_phenotypes/Transcriptomics/Discovery_cohort_RNAseq")

# Output directories
output_dir_2 <- "Outputs/Files/ORA"
dir.create(output_dir_2, recursive = TRUE, showWarnings = FALSE)

# Load cluster DEA results
cluster_files <- list.files(
  "./Outputs/Files/Limma",
  pattern = "Cluster_\\d+_DEGs_all.csv",
  full.names = TRUE
)

# Create a named list of DEA results
dea_results_list <- lapply(cluster_files, read.csv, row.names = 1)
names(dea_results_list) <- paste0("Cluster", seq_along(dea_results_list))

# Parameters for DEG selection
padj_cutoff <- 0.01
logfc_cutoff <- 0   # set to e.g. 1 if you want directionality

# Background genes (union of all tested genes)
background_genes <- unique(
  unlist(lapply(dea_results_list, rownames))
)

background_genes <- toupper(background_genes)

# Gene list per cluster
gene_list <- lapply(dea_results_list, function(df) {
  
  df <- as.data.frame(df)
  
  sig_genes <- df %>%
    rownames_to_column("gene") %>%
    filter(
      adj.P.Val < padj_cutoff,
      logFC > logfc_cutoff
    ) %>%
    arrange(adj.P.Val) %>%        # most significant first
    slice_head(n = 100) %>%       # keep at most top 100
    pull(gene)
  
  toupper(sig_genes)
})

# Name the list consistently
names(gene_list) <- names(dea_results_list)

# Quick sanity check
sapply(gene_list, length)

# MSigDB C2 gene sets
msigdb_c2 <- msigdbr(species = "Homo sapiens", category = "C2")
msigdb_c2_cp <- msigdb_c2 %>%
  filter(gs_subcollection %in% c("CP:KEGG","CP:REACTOME","CP:BIOCARTA",
                                 "CP:WIKIPATHWAYS", "CP:PID", "CGP:NABA"))

gene_sets <- msigdb_c2_cp[, c("gs_name", "gene_symbol")]
gene_sets$gene_symbol <- toupper(gene_sets$gene_symbol)

## load the pathway analysis function
source("02_analysis_and_figure_scripts/pathway_analysis_function.R")

## perform pathway analysis for each cluster of genes
for (cluster in names(gene_list)) {
  pathway_analysis(
    gene_list        = gene_list[[cluster]],
    background_genes = background_genes,
    gene_sets        = gene_sets,
    csv_file         = sprintf("Outputs/Files/ORA/rna_%s_enriched_pathways.csv", cluster)
  )
}
