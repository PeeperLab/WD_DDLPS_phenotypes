

set.seed(123)

## ---- Working directory ----------------------------------------------------
## Before running, set the working directory to the folder holding the
## raw RNA-seq input data (expected folder layout: see README), e.g.
## setwd("/path/to/RNAseq_raw_data")
library(readxl)
library(dplyr)

meta_file <- "meta/raw_metadata/aspera/2025-10-30-metadata-cp1240-2025102847001081-for-Customer.tsv"
file2 <- "meta/raw_metadata/aspera/DKTK-P157_Clinical_Information_251202_TJ3.xlsx"

out_dir1 <- "meta/metadata"
if (!dir.exists(out_dir1)) {
  dir.create(out_dir1, recursive = TRUE)
}

raw_metadata <- read.delim(meta_file)
raw_metadata <- as.data.frame(raw_metadata)
dim(raw_metadata)
length(unique(raw_metadata$SAMPLE_NAME))

table(raw_metadata$SAMPLE_TYPE)
rna_meta <- raw_metadata[raw_metadata$SEQUENCING_TYPE == "RNA", ]
rna_meta$sample_id <- sub("-R[12].*$", "", rna_meta$SAMPLE_NAME)
rna_meta$fastq_id <- sub("^(([^-]*-){1}[^-]*).*", "\\1", rna_meta$FASTQ_FILE)
rna_meta$sample_type <- gsub("[0-9]", "", rna_meta$SAMPLE_TYPE)

samples <- unique(rna_meta$SAMPLE_NAME)

sub_rna_meta <- rna_meta[,c("sample_id","fastq_id","sample_type","INSTRUMENT_MODEL",
                            "SAMPLE_TYPE","TISSUE_TYPE","LIB_PREP_KIT", "SEQUENCING_KIT",
                            "BASE_MATERIAL","CENTER_NAME","SEQUENCING_TYPE","INSTRUMENT_PLATFORM")]

sub_rna_meta <- unique(sub_rna_meta)

# write.csv(sub_rna_meta,
#           file.path(out_dir1, "selected_metadata.csv"),
#           row.names = FALSE)

metadata1 <- sub_rna_meta
# patient = between the first and second "-"
metadata1$patient <- sub("^H021-([^-]+)-.*$", "\\1", metadata1$sample_id)
# move patient to the first column
metadata1 <- metadata1[, c("patient", setdiff(names(metadata1), "patient"))]

metadata2 <- read_excel(file2)
metadata2 <- as.data.frame(metadata2)
metadata2 <- metadata2[metadata2$`Sequencing data` == "yes",]
metadata2 <- metadata2[metadata2$`Sequencing data` == "yes",]
metadata2[metadata2$tumor_type_TJ2 == "LPS", "tumor_type_TJ2"] <- "WD/DDLPS"
metadata2[metadata2$tumor_type_TJ2 == "WD+DDLPS", "tumor_type_TJ2"] <- "WD/DDLPS"


sub_meta2 <- metadata2[,c("MASTER pid","tumor_type_TJ2")]


metadata1 <- metadata1 %>%
  left_join(
    metadata2 %>%
      transmute(patient = as.character(`MASTER pid`),
                tumor_type = as.character(tumor_type_TJ2),
                gender = as.character(Sex),
                Deathdate = Deathdate,
                Deathdate_accuracy = `Deathdate Accuracy`,
                Date_diagnosis = `Date of Diagnosis`,
                last_contact_date = `Date of last contact`) %>%
      distinct(patient, .keep_all = TRUE),
    by = "patient"
  )

metadata1$dataset <- rep("Aspera",nrow(metadata1))
metadata1[metadata1$fastq_id == "AS-256974", "sample_type"] <- "Fat_tissue"
metadata1[metadata1$fastq_id == "AS-256974", "tumor_type"] <- "Fat_tissue"
# 
# write.csv(metadata1,
#           file.path("combined_rna_metadata.csv"),
#           row.names = FALSE)
# 
# 
# head(metadata1)
sub_metadata <- metadata1[,c("patient","sample_id","fastq_id","tumor_type","dataset","INSTRUMENT_MODEL","sample_type","gender","Deathdate","Date_diagnosis","last_contact_date","Deathdate_accuracy")]

colnames(sub_metadata) <- c("petient","sampleID","fastq_id", "subtype","dataset","seq_tech","tissue_type","gender","Deathdate","Date_diagnosis","last_contact_date","Deathdate_accuracy")


write.csv(sub_metadata,
          file.path(out_dir1,"aspera_final_rna_metadata.csv"),
          row.names = FALSE)


#### calculate surival days and add status column
meta2 <- sub_metadata %>%
  mutate(
    across(
      c(Date_diagnosis, Deathdate, last_contact_date),
      ~ as.Date(na_if(., "NA"), format = "%d.%m.%Y")
    ),
    
    ##calculate survival days
    surv_time_days = as.numeric(
      if_else(
        !is.na(Deathdate),
        Deathdate - Date_diagnosis,        
        last_contact_date - Date_diagnosis 
      )
    ),
    
    ## add status column
    status = if_else(!is.na(Deathdate), 1L, 0L)
  )

write.csv(meta2,
          file.path(out_dir1,"aspera_rna_metadata_with_survival.csv"),
          row.names = FALSE)