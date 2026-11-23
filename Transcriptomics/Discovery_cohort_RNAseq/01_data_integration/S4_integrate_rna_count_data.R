###This script is used to integrate the rna-seq counts data.
library(readr)
library(tidyr)
library(tibble)
library(stringr)
library(dplyr)
library(readxl)

# set working directory
## ---- Working directory ----------------------------------------------------
## Before running, set the working directory to the folder holding the
## raw RNA-seq input data (expected folder layout: see README), e.g.
## setwd("/path/to/RNAseq_raw_data")
## file path
file_path1 <- "meta/metadata/sarc_lps_metadata.csv"
file_path2 <- "data/rds/Song2022.xlsx"
file_path3 <- "data/rds/bevil2023_count_matrix.txt"
file_path4 <- "data/rds/GSE221492_counts_matrix.tsv"
file_path5 <- "data/rds/sarc_lps_count_data.rds"
file_path6 <- "meta/metadata/metadata_282597.csv"
file_path7 <- "meta/metadata/metadata_987378.csv"
file_path8 <- "data/rds/counts_data_282597.rds"
file_path9 <- "data/rds/counts_data_987378.rds"
file_path10 <-"meta/metadata/metadata_1053633.csv"
file_path11 <- "data/rds/counts_data_1053633.rds"
file_path12 <- "meta/metadata/aspera_final_rna_metadata.csv"
file_path13 <- "data/rds/counts_data_german.rds"

out_dir1   <- "meta/metadata" 
out_dir2 <- "data/rds"

### 1. integrate metadata
#TCGA_SARC metadata
sarc_metadata <- read_csv(file_path1)
sarc_metadata$seq_tech <- rep("Illumina HiSeq 2000",nrow(sarc_metadata))
dim(sarc_metadata)
## song dataset
# tumor sheet
song_tumor_data <- as.data.frame(read_excel(file_path2, sheet = "Supplementary2", col_names = FALSE))
tumor_metadata <- song_tumor_data[c(3,4),]
song_metadata_tumor <- as.data.frame(t(tumor_metadata))
colnames(song_metadata_tumor) <- song_metadata_tumor[1,]
song_metadata_tumor <- song_metadata_tumor[-1,]
colnames(song_metadata_tumor) <- c("sampleID", "subtype")

#control sheet
song_control_data <- read_excel(file_path2, sheet = "Supplementary3", col_names = FALSE)
control_metadata <- as.data.frame(song_control_data[3,])
rownames(control_metadata) <- "sampleID"
new_row <- rep("normal_control", ncol(control_metadata))
control_metadata <- rbind(control_metadata[1, ], new_row,control_metadata[-1, ])
rownames(control_metadata)[2] <- "subtype"

song_metadata_control <-  as.data.frame(t(control_metadata))
song_metadata_control <- song_metadata_control[-1,]

# combine the tumor and contol metadata
song_metadata <- rbind(song_metadata_tumor,song_metadata_control)
song_metadata$dataset <- rep("Song2022",nrow(song_metadata))
song_metadata$seq_tech <- rep("Illumina ",nrow(song_metadata))
rownames(song_metadata) <- NULL
head(song_metadata)
table(song_metadata$subtype)
dim(song_metadata)
##bevil data
bevil_data <- read.table(file_path3)

bevil_metadata <- data.frame(sampleID=colnames(bevil_data))
bevil_metadata$subtype <- ifelse(grepl("DD", bevil_metadata$sampleID), "DDLPS",
                                 ifelse(grepl("WD", bevil_metadata$sampleID), "WDLPS",
                                        ifelse(grepl("N", bevil_metadata$sampleID), "Fat_tissue", NA)))

bevil_metadata$dataset <- rep("Bevill2023",nrow(bevil_metadata))
bevil_metadata$seq_tech <- rep("Illumina NextSeq 500",nrow(bevil_metadata))
head(bevil_metadata)
## gse221492
gse221492 <- read_tsv(file_path4)

gse221492_clean <- as.data.frame(gse221492)

# delete the columns with "Hibernoma" or "lipoma" in their colunm names
gse221492_clean <- gse221492_clean %>%
  dplyr::select(-matches("Hibernoma|lipoma"))

gse221492_metadata <- data.frame(sampleID=colnames(gse221492_clean))
gse221492_metadata <- as.data.frame(gse221492_metadata[-c(1,2),])
colnames(gse221492_metadata) <- "sampleID"
# extract the content after the last "_" in sampleID 
sample_last_part <- str_extract(gse221492_metadata$sampleID, "[^_]+$")

# create the subtype based on the extracted content
subtype <- #ifelse(grepl("WDLPS",sample_last_part),"WDLPS",
  ifelse(grepl("DD",sample_last_part),"DDLPS",
         ifelse(grepl("WD",sample_last_part),"WDLPS",
                ifelse(grepl("fat",sample_last_part),"Fat_tissue", NA)))

gse221492_metadata$subtype <- subtype

gse221492_metadata$dataset <- rep("GSE221492",nrow(gse221492_metadata))
gse221492_metadata$seq_tech <- rep("Illumina NovaSeq 6000",nrow(gse221492_metadata))

## 282597
metadata_282597 <- read_csv(file_path6)
metadata_282597 <- as.data.frame(metadata_282597)
metadata_282597$seq_tech <- rep("Illumina HiSeq 2000",nrow(metadata_282597))

##987378
metadata_987378 <- read_csv(file_path7)
metadata_987378 <- as.data.frame(metadata_987378)
metadata_987378$seq_tech <- rep("Illumina NovaSeq 6000",nrow(metadata_987378))

## 53633
metadata_53633 <- read_csv(file_path10)
metadata_53633 <- as.data.frame(metadata_53633)
metadata_53633$seq_tech <- rep("Illumina HiSeq 4000",nrow(metadata_53633))
### asprea
metadata_aspera <- read.csv(file_path12)
metadata_aspera <- as.data.frame(metadata_aspera)
metadata_aspera <- metadata_aspera[ ,colnames(metadata_aspera) %in% c("fastq_id","subtype", "dataset", "seq_tech")]

colnames(metadata_aspera) <- c("sampleID","subtype", "dataset", "seq_tech")

### combine th emetadata from 8 datasets
integrate_metadata <- rbind(sarc_metadata,song_metadata,bevil_metadata,gse221492_metadata,metadata_282597,metadata_987378, metadata_53633,metadata_aspera)
integrate_metadata <- integrate_metadata %>%
  mutate(subtype = case_when(subtype == "pleomorphic" ~ "PLPS",
                           subtype == "normal_control" ~ "Fat_tissue",
                           TRUE ~ subtype))

#integrate_metadata$subtype <- factor(integrate_metadata$subtype, levels = c("DDLPS", "WDLPS", "Fat_tissue", "PLPS"))
#integrate_metadata <- integrate_metadata[order(integrate_metadata$subtype), ]
dim(integrate_metadata)
#integrate_metadata <- integrate_metadata[!integrate_metadata$subtype == "PLPS", ]
#dim(integrate_metadata)
# Dataset labels used in the manuscript (accession numbers)
dataset_labels <- c(Aspera     = "EGAS00001004813",
                    GSE221492  = "PRJNA914687",
                    Bevill2023 = "PRJNA880357",
                    Song2022   = "Song et al")
idx <- integrate_metadata$dataset %in% names(dataset_labels)
integrate_metadata$dataset[idx] <- unname(dataset_labels[integrate_metadata$dataset[idx]])
table(integrate_metadata$dataset)

write_csv(integrate_metadata,file = file.path(out_dir1,"Metadata_discovery_all.csv"))


###The following part is used to integrate gene counts data.
sarc_count_data <- readRDS(file_path5)
# check if there are rows with row sum =0
num_zero_expression_rows <- sum(rowSums(sarc_count_data[, 2:ncol(sarc_count_data)]) == 0)

# delecte the rows with rowsum=0
nonsum0_count_data <- sarc_count_data[rowSums(sarc_count_data[, 2:ncol(sarc_count_data)]) != 0, ]

#for duplicate genes, use the sum of each duplicate genes
sarc_count_data_cleaned <- nonsum0_count_data %>%
  group_by(gene_symbol) %>%
  summarise(across(starts_with("TCGA"), ~ sum(.x, na.rm = TRUE)))

### song data
# gene counts data of tumor samples
song_tumor_count <- song_tumor_data[-c(1,2,4),]
colnames(song_tumor_count) <- song_tumor_count[1,]
song_tumor_count <- song_tumor_count[-1,]
# calculate the number of rowSum = 0
song_tumor_count[, 2:ncol(song_tumor_count)] <- lapply(song_tumor_count[, 2:ncol(song_tumor_count)], as.numeric)
num_sum0 <- sum(rowSums(song_tumor_count[,2:ncol(song_tumor_count)])==0)
print(num_sum0)

# gene counts data of control samples
song_control_count <- song_control_data[-c(1,2),]
colnames(song_control_count) <- song_control_count[1, ]
colnames(song_control_count)[is.na(colnames(song_control_count))] <- "Gene.Symbol"
song_control_count <- song_control_count[-1,]
# calculate the number of rowSum = 0
song_control_count[, 2:ncol(song_control_count)] <- lapply(song_control_count[, 2:ncol(song_control_count)], as.numeric)
num_sum0 <- sum(rowSums(song_control_count[,2:ncol(song_control_count)])==0)
print(num_sum0)

# combine the tumor counts data with the control data
song_count_data <- inner_join(song_tumor_count, song_control_count, by = "Gene.Symbol")
# seperate the ensembl id and gene symbol
song_count_data <- song_count_data %>% 
  separate(Gene.Symbol, into = c("gene_id", "gene_symbol"), sep = "\\|")
song_count_data <- song_count_data[,-1]
head(song_count_data)
dim(song_count_data)

# use the sum of each duplicate genes
song_count_data_cleaned <- song_count_data %>%
  group_by(gene_symbol) %>%
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)))

## merge TCGA and song dataset
merged_count_data1 <- inner_join(sarc_count_data_cleaned,song_count_data_cleaned, by = "gene_symbol" )

### bevil data
bevil_count_data <- bevil_data %>%
  rownames_to_column(var="gene_symbol")
#check if there are empty gene_symbol
print(sum(bevil_count_data$gene_symbol==""))
# check if there are NAs in the dataset
any(is.na(bevil_count_data))
# only keep the rows with rowSum is not 0
bevil_count_data_nonsum0 <- bevil_count_data[rowSums(bevil_count_data[,2:ncol(bevil_count_data)])!=0,]

#check if there are duplicated genes
gene_counts3 <- table(bevil_count_data_nonsum0$gene_symbol)
# filter the duplicate genes
duplicate_genes3 <- gene_counts3[gene_counts3 > 1]
print(length(duplicate_genes3))

## merge bevil data and merged data1
merged_count_data2 <- inner_join(merged_count_data1, bevil_count_data_nonsum0, by = "gene_symbol" )

####gse221492
gse221492_count_data <- gse221492_clean[,-1]

colnames(gse221492_count_data)[colnames(gse221492_count_data)=="hgnc_symbol"] <- "gene_symbol"
#check if there are empty gene_symbol
print(sum(gse221492_count_data$gene_symbol==""))
# check if there are NAs in the dataset
any(is.na(gse221492_count_data))
#delete the rows with NAs
gse221492_count_data <- na.omit(gse221492_count_data)
# check if there are still NAs in the dataset
any(is.na(gse221492_count_data))
# only keep the rows with rowSum is not 0
gse221492_count_data <- gse221492_count_data[rowSums(gse221492_count_data[,2:ncol(gse221492_count_data)])!=0,]

# use the sum of each duplicate genes
gse221492_count_data_cleaned <- gse221492_count_data %>%
  group_by(gene_symbol) %>%
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)))

### merge gse221492 data and merged data2
merged_count_data3 <- inner_join(merged_count_data2, gse221492_count_data_cleaned, by = "gene_symbol" )


###282597
count_data_282597 <- readRDS(file_path8)
### merge 282597 data and merged data3
merged_count_data4 <- inner_join(merged_count_data3, count_data_282597, by = "gene_symbol" )

###987378
count_data_987378 <- readRDS(file_path9)

### merge 987378 data and merged data4
merged_count_data5 <- inner_join(merged_count_data4, count_data_987378, by = "gene_symbol" )

###53633
count_data_53633 <- readRDS(file_path11)

### merge 53633 data and merged data5
merged_count_data6 <- inner_join(merged_count_data5, count_data_53633, by = "gene_symbol" )
dim(merged_count_data6)

### asprea dataset
count_aspera <- readRDS(file_path13)
head(count_aspera)
any(is.na(count_aspera))
dim(count_aspera)
### merge asprea dataset and merged data6
merged_count_data7 <- inner_join(merged_count_data6, count_aspera, by = "gene_symbol" )
dim(merged_count_data7)
head(merged_count_data7)
merged_count_data7 <- as.data.frame(merged_count_data7)
merged_count_data7 <- merged_count_data7[,c("gene_symbol",integrate_metadata$sampleID)]
dim(merged_count_data7)
table(integrate_metadata$subtype)
### save the merged count data
saveRDS(merged_count_data7,file = file.path(out_dir2,"integrated_rna_count_data.rds") )
