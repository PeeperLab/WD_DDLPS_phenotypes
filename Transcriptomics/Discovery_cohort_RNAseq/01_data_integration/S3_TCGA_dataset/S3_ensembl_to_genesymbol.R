### this script is used to convert the ensembl id to gene symbol

set.seed(123)
## ---- Working directory ----------------------------------------------------
## Before running, set the working directory to the folder holding the
## raw RNA-seq input data (expected folder layout: see README), e.g.
## setwd("/path/to/RNAseq_raw_data")

library(biomaRt)
library(dplyr)
library(tibble)

input <- "data/rds/sarc_lps_count_ensembl.rds"
data <- readRDS(input)
data_t <- as.data.frame(t(data))
colnames(data_t) <- data_t[1,]
cleaned_data_t <- data_t[-1,]
cleaned_data_t <- cleaned_data_t %>%
  rownames_to_column(var = "gene_id")
#delete the version number ".*"
cleaned_data_t$gene_id <- gsub("\\..*", "", cleaned_data_t$gene_id)
#gene ids should be character
gene_ids <- as.character(cleaned_data_t[[1]])  
#connect Ensembl dataset
ensembl <- useEnsembl(biomart = "ensembl", 
                      dataset = "hsapiens_gene_ensembl")

gene_info <- getBM(
  attributes = c("ensembl_gene_id","hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = gene_ids,
  mart = ensembl
)

#rename the colnames
colnames(gene_info) <- c("gene_id","gene_symbol")
sarc_raw_count_data <- inner_join(gene_info, cleaned_data_t, by = "gene_id")

# delete the rows without gene_symbol
sarc_nonNA_count_data <- sarc_raw_count_data[sarc_raw_count_data$gene_symbol!="",]
dim(sarc_nonNA_count_data)

# calculate the number of rowSum = 0
sarc_nonNA_count_data[, 3:ncol(sarc_nonNA_count_data)] <- lapply(sarc_nonNA_count_data[, 3:ncol(sarc_nonNA_count_data)], as.numeric)

num_zero_expression_rows <- sum(rowSums(sarc_nonNA_count_data[, 3:ncol(sarc_nonNA_count_data)]) == 0)
print(num_zero_expression_rows)
# kepp the rows that tha rowsum is not 0
sarc_non0_count_data <- sarc_nonNA_count_data[rowSums(sarc_nonNA_count_data[, 3:ncol(sarc_nonNA_count_data)]) != 0, ]

sarc_non0_count_data <- sarc_non0_count_data[,-1]
head(sarc_non0_count_data)

saveRDS(sarc_non0_count_data,"data/rds/sarc_lps_count_data.rds")
