set.seed(123)
## ---- Working directory ----------------------------------------------------
## Before running, set the working directory to the folder holding the
## raw WGS input data (expected folder layout: see README), e.g.
## setwd("/path/to/WGS_raw_data")

library(readxl)
library(dplyr)
library(data.table)
meta_file <- "WGS_input/2025-10-30-metadata-cp1240-2025102847001081-for-Customer.tsv"
meta_file2 <- "WGS_input/DKTK-P157_Clinical_Information_251202_TJ3.xlsx"
cnv_dir   <- "WGS_input/renamed/common/cnv"
out_dir <- "WGS_input/metadata"

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

##### handle metadata

raw_metadata <- read.delim(meta_file)
raw_metadata <- as.data.frame(raw_metadata)

dna_meta <- raw_metadata[raw_metadata$SEQUENCING_TYPE != "RNA", ]
dna_meta$sample_id <- sub("-D[12].*$", "", dna_meta$SAMPLE_NAME)
dna_meta$patient <- sub("^H021-([^-]+)-.*$", "\\1", dna_meta$sample_id)
dna_meta$fastq_id <- sub("^(([^-]*-){1}[^-]*).*", "\\1", dna_meta$FASTQ_FILE)

dna_meta$sample_type <- gsub("[0-9]", "", dna_meta$SAMPLE_TYPE)
table(dna_meta$sample_type)
# write.csv(dna_meta,
#           file.path(out_dir, "raw_dna_metadata.csv"),
#           row.names = FALSE)


samples <- unique(dna_meta$SAMPLE_NAME)


sub_dna_meta <- dna_meta[,c("patient","sample_id","fastq_id","sample_type","SEQUENCING_TYPE","INSTRUMENT_MODEL",
                            "SAMPLE_TYPE","TISSUE_TYPE","LIB_PREP_KIT", "SEQUENCING_KIT",
                            "BASE_MATERIAL","CENTER_NAME","INSTRUMENT_PLATFORM")]

sub_dna_meta <- unique(sub_dna_meta)

# write.csv(sub_dna_meta,
#           file.path(out_dir1, "selected_dna_metadata.csv"),
#           row.names = FALSE)



### clinical information
metadata2 <- read_excel(meta_file2 ,sheet = "Clinical Information")
metadata2 <- as.data.frame(metadata2)
metadata2 <- metadata2[metadata2$`Sequencing data` == "yes",]
metadata2[metadata2$tumor_type_TJ2 == "LPS", "tumor_type_TJ2"] <- "WD/DDLPS"
metadata2[metadata2$tumor_type_TJ2 == "WD+DDLPS", "tumor_type_TJ2"] <- "WD/DDLPS"

select_dna_meta <- sub_dna_meta[sub_dna_meta$sample_type %in% c("metastasis","tumor","undefined-neoplasia"),]
select_dna_meta[select_dna_meta$sample_type =="undefined-neoplasia","sample_type"] <- "tumor"
select_dna_meta <- select_dna_meta[!duplicated(select_dna_meta$sample_id), ]

dna_meta <- select_dna_meta %>%
  left_join(
    metadata2 %>%
      transmute(patient = as.character(`MASTER pid`),
                tumor_type = as.character(tumor_type_TJ2),
                gender = as.character(Sex),
                Deathdate = Deathdate,
                Deathdate_accuracy = `Deathdate Accuracy`,
                Date_diagnosis = `Date of Diagnosis`,
                last_contact_date = `Date of last contact`
                ) %>%
      distinct(patient, .keep_all = TRUE),
    by = "patient"
  )


write.csv(dna_meta,
          file.path(out_dir, "combined_dna_metadata_with_survival.csv"),
          row.names = FALSE)


#### extract metadata from file names

files <- list.files(
  cnv_dir,
  pattern = ".(txt|tsv)$",
  full.names = TRUE,
  recursive = TRUE
)
if (length(files) == 0) stop("No files found in ", cnv_dir)


parse_aceseq_header <- function(f) {
  hdr <- readLines(f, n = 50, warn = FALSE)
  get_val <- function(prefix) {
    line <- hdr[grepl(paste0("^#", prefix, ":"), hdr)]
    if (length(line) == 0) return(NA_real_)
    as.numeric(sub(paste0("^#", prefix, ":"), "", line[1]))
  }
  data.table(
    sample  = sub("\\.(txt|tsv)$", "", basename(f), ignore.case = TRUE),
    ploidy  = get_val("ploidy"),
    tcc     = get_val("tcc"),
    quality = get_val("quality")
  )
}


meta <- rbindlist(lapply(files, parse_aceseq_header), fill = TRUE)
meta[, sample := make.unique(sample)]
meta[, patient := sub("^[^-]+-([^-]+).*", "\\1", sample)]
setcolorder(meta, c("patient", setdiff(names(meta), "patient")))
meta <- as.data.frame(meta)

## combine subtype information of each patient
dna_subtype <- dna_meta %>%
  dplyr::mutate(
    patient    = as.character(patient),
    tumor_type = as.character(tumor_type),
    tissue_type = as.character(sample_type),
    gender = as.character(gender)
  ) %>%
  dplyr::distinct(patient, .keep_all = TRUE) %>%
  dplyr::select(patient, tumor_type,tissue_type,gender)

cna_meta <- as.data.frame(meta) %>%
  dplyr::mutate(patient = as.character(patient)) %>%
  dplyr::left_join(dna_subtype, by = "patient") 


write.csv(cna_meta, file.path(out_dir, "final_dna_metadata.csv"), row.names = FALSE)


