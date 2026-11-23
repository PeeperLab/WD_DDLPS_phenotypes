
###This script is used to handle TCGA-SARC dataset and select the survival data.

# set working directory
set.seed(123)
## ---- Working directory ----------------------------------------------------
## Before running, set the working directory to the folder holding the
## raw RNA-seq input data (expected folder layout: see README), e.g.
## setwd("/path/to/RNAseq_raw_data")

library(readr)
library(dplyr)
library(tibble)


# clinical data path
file_path1 <- "data/TCGA/sarc_clin_data.csv"
# raw count data path
file_path2 <- "data/TCGA/sarc_raw_data.csv"

out_dir1   <- "meta/metadata" 
out_dir2 <- "data/rds"

# handle the clinical data
clin_data <- read_csv(file_path1)
# select liposarcoma data
lps_data <- subset(clin_data, grepl("liposarcoma", primary_pathology_histological_type, ignore.case = TRUE))

# extract the bcr_patient_barcode and primary_pathology_histological_type columns
filtered_lps_data <- lps_data[, c("bcr_patient_barcode", "primary_pathology_histological_type")]
colnames(filtered_lps_data) <- c("sampleID","subtype")


# process the count data
count_data <- read_csv(file_path2)
# transpose the count data to make sure the sample ids are the row names
count_data_t <- as.data.frame(t(count_data))
colnames(count_data_t) <- count_data_t[1,]
count_data_t <- count_data_t[-1,] %>%
  rownames_to_column(var="sampleid")

#create a new column to store the first 12 characters of sample id
count_data_t <- count_data_t %>%
  mutate(sampleid_prefix = substr(sampleid, 1, 12)) %>%
  dplyr::select(sampleid_prefix, everything())


# check if tehre are samples only in filtered_lps_data1.
diff1 <- setdiff(filtered_lps_data$sampleID, count_data_t$sampleid_prefix)

filtered_lps_data <- filtered_lps_data[filtered_lps_data$sampleID != diff1, ]


# combine these two datasets by  bcr_sample_barcode == sampleid_prefix 
merged_data <- left_join(filtered_lps_data,count_data_t, by = c("sampleID"="sampleid_prefix"))

# check if there are duplicated sampleIDs
duplicated_rows <- merged_data %>%
  filter(duplicated(sampleID) | duplicated(sampleID, fromLast = TRUE))

#TCGA-FX-A2QS-11A is normal tissue, delete this row
#TCGA-FX-A2QS-01A is tumor, keep this row
cleaned_data <- merged_data %>%
  filter(!grepl("TCGA-FX-A2QS-11A", sampleid))

# select sample Id and subtype columns to generate the lps metadata
sarc_metadata <- cleaned_data[,c(1,2)]
sarc_metadata$subtype <- "DDLPS"
sarc_metadata$dataset <- rep("TCGA_SARC",nrow(sarc_metadata))

#save the sarc lps metadata
write_csv(sarc_metadata,file.path(out_dir1,"sarc_lps_metadata.csv"))

### select the survival data
survival_data <- lps_data [lps_data$bcr_patient_barcode %in% cleaned_data$sampleID, ]
survival_data <- survival_data[,c("bcr_patient_barcode", "days_to_death", "days_to_last_followup")]
colnames(survival_data)[1] <-  "sampleID"
write_csv(survival_data,file.path(out_dir1,"TCGA_sarc_dd_survival_data.csv"))


# only keep sample id and count data
sarc_lps_count_data <- cleaned_data[,-c(2,3)]

# save the count data, this file will be used for following analysis
saveRDS(sarc_lps_count_data,file.path(out_dir2,"sarc_lps_count_ensembl.rds"))






