### This script uses limma to normalize RNA-seq count data and combat to correct batch effect.
set.seed(123)
# set working directory
## ---- Working directory ----------------------------------------------------
## Before running, set the working directory to the folder holding the
## raw RNA-seq input data (expected folder layout: see README), e.g.
## setwd("/path/to/RNAseq_raw_data")
library(readr)
library(sva)
library(limma)
library(edgeR)

# input data path
file_path1 <-"meta/metadata/Metadata_discovery_all.csv"
file_path2 <- "data/rds/integrated_rna_count_data.rds"

metadata <- read_csv(file_path1)
count_data <- readRDS(file_path2)

count_data_unique <- as.data.frame(count_data)
# change the gene_symbol column to row names
rownames(count_data_unique) <- count_data_unique$gene_symbol
count_data_unique$gene_symbol <- NULL
setdiff(metadata$sampleID,colnames(count_data_unique) )
count_data_unique <- count_data_unique[,metadata$sampleID]
all(metadata$sampleID == colnames(count_data_unique))

count_data <- count_data_unique

### data normalization
exprSet <- count_data
metadata$subtype <- factor(metadata$subtype)
subtype <- metadata$subtype

ncol(exprSet) == length(subtype)
# wrap the raw count matrix and sample groups into an edgeR object
dge <- DGEList(counts = exprSet, group = subtype)

# Calculate counts per million (CPM) for each gene in dge
cpm_values <- cpm(dge)

# filter out low expression genes
keep <- filterByExpr(dge, group = subtype)   
dge <- dge[keep, , keep.lib.sizes = FALSE]
#perform TMM normalization, correcting for differences in library size and composition
dge <- calcNormFactors(dge, method = "TMM")
# buld design matrix: rows are samples, columns are subtypes, 1 means this sample belongs to this subtype
design <- model.matrix(~ 0 + subtype,data=metadata)
colnames(design) <- levels(subtype)
rownames(design) <- colnames(exprSet)

# transform the normalized counts into log₂‑CPM values, model the mean–variance relationship, and generate precision weights 
voom_data <- voom(dge, design, plot = T)
# get the normalized log2 cpm data
normalized_data <- as.data.frame(voom_data$E)

#### remove batch effect
metadata$seq_tech <- factor(metadata$seq_tech)
metadata$dataset  <- factor(metadata$dataset)
metadata$subtype  <- factor(metadata$subtype)

batch_joint <- interaction(metadata$seq_tech, metadata$dataset, drop = TRUE)
mod <- model.matrix(~ subtype, data = metadata)

combat_data <- ComBat(
  dat  = as.matrix(normalized_data),
  batch = batch_joint,
  mod  = mod,
  par.prior = TRUE,
  prior.plots = FALSE
)


## save the batch corrected data
combat_data_df <- as.data.frame(combat_data)
dim(combat_data_df)
saveRDS(combat_data_df,"data/rds/Batch_corrected_data_discovery_all.rds")

