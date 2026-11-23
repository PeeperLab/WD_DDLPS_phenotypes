###This script is used to download TCGA DATA.
set.seed(123)
### Set working directory
## ---- Working directory ----------------------------------------------------
## Before running, set the working directory to the folder holding the
## raw RNA-seq input data (expected folder layout: see README), e.g.
## setwd("/path/to/RNAseq_raw_data")

library(TCGAbiolinks)
library(SummarizedExperiment)

projects <- getGDCprojects()

download_directory <- "data/TCGA"
dir.create(download_directory)

# download the count data
query_exp <- GDCquery(project = 'TCGA-SARC',
                      data.category = 'Transcriptome Profiling',
                      data.type = "Gene Expression Quantification",
                      experimental.strategy = 'RNA-Seq'
)
output_query <- getResults(query_exp)

# download the data
GDCdownload(query_exp, directory = download_directory)
tcga_srac <- GDCprepare(query_exp, directory = download_directory, summarizedExperiment = TRUE)

# get the raw count data
sarc_matrix <- assay(tcga_srac,"unstranded")

# Save the matrix to a CSV file
write.csv(sarc_matrix, file = file.path(download_directory,"sarc_raw_data.csv"), row.names = TRUE)


# download the clinical data
query_clin <- GDCquery(project = "TCGA-SARC",
                       data.category = "Clinical",
                       data.type = "Clinical Supplement",
                       data.format = "BCR XML")
output_clin <- getResults(query_clin)

# download the data
GDCdownload(query_clin)

# prepare the clinical data
clinical_data <- GDCprepare_clinic(query_clin, clinical.info = "patient")

write.csv(clinical_data, file = file.path(download_directory,"sarc_clin_data.csv"), row.names = TRUE)
