# This script performs ssGSEA on the proteomic dataset using the GSVA package.

library(readr)
library(dplyr)
library(tibble)
library(msigdbr)
library(GSVA)

set.seed(1507)
## ---- Working directory ----------------------------------------------------
## Before running, set the working directory to the "Proteomics" folder of
## this repository (the folder that holds Inputs/ and Outputs/), e.g.
## setwd("/path/to/WD_DDLPS_phenotypes/Proteomics")

output_dir <- "Outputs/Files/ssGSEA"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

## 1. prepare the data
# expression data
exprs_data <- read.csv("Inputs/Proteomic_data_batch_corrected.csv", row.names = 1)
exprs_data <- t(exprs_data)

#RNA signatures:
pathway_list <- c(
  "PID_TCR_PATHWAY",
  "WP_MODULATORS_OF_TCR_SIGNALING_AND_T_CELL_ACTIVATION",
  "WP_B_CELL_RECEPTOR_SIGNALING",
  "NAKAJIMA_MAST_CELL",
  "REACTOME_NEUTROPHIL_DEGRANULATION",
  "REACTOME_ANTIGEN_PROCESSING_CROSS_PRESENTATION",
  "REACTOME_MHC_CLASS_II_ANTIGEN_PRESENTATION",
  "REACTOME_G2_M_CHECKPOINTS",
  "REACTOME_CELL_CYCLE_CHECKPOINTS",
  "REACTOME_MITOTIC_SPINDLE_CHECKPOINT",
  "REACTOME_COLLAGEN_CHAIN_TRIMERIZATION",
  "REACTOME_COLLAGEN_DEGRADATION",
  "REACTOME_COLLAGEN_FORMATION",
  "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION",
  "REACTOME_NCAM1_INTERACTIONS",
  "WP_NEURAL_CREST_DIFFERENTIATION",
  "REACTOME_ADIPOGENESIS",
  "REACTOME_FATTY_ACID_METABOLISM",
  "WP_PPAR_SIGNALING",
  "WP_DIFFERENTIATION_OF_WHITE_AND_BROWN_ADIPOCYTE"
)

# load ALL C2 canonical pathways (KEGG, Reactome, NABA, PID, etc.)
msigdb_c2_sets <- msigdbr(species = "Homo sapiens", category = NULL, subcategory = NULL)

# Filter for your manually defined pathways
filtered_gene_sets <- msigdb_c2_sets %>%
  filter(gs_name %in% pathway_list)

# convert to gene list: the input for ssGSEA
ssgsea_gene_list <- split(filtered_gene_sets$gene_symbol, filtered_gene_sets$gs_name)
#print(ssgsea_gene_list)

## single sample GSEA
exprs_matrix <- as.matrix(exprs_data)

ssgsea_par <- ssgseaParam(exprs_matrix, 
                          ssgsea_gene_list,
                          minSize = 5, 
                          maxSize = 500,
                          normalize = TRUE)

gsva_results <- gsva(ssgsea_par)
# Keep only pathways that exist in the gsva_results rownames
available_pathways <- intersect(pathway_list, rownames(gsva_results))

# Subset safely
gsva_results <- gsva_results[available_pathways, ]
gsva_results <- as.data.frame(gsva_results) %>%
  rownames_to_column(var = "pathway")

# Save CSV
write.csv(
  gsva_results,
  file = file.path(output_dir, "gsva_results_RNA_sigs.csv"),
)

