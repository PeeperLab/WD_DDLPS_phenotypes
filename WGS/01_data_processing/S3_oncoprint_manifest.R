
set.seed(123)

### Working directory and inputs
## ---- Working directory ----------------------------------------------------
## Before running, set the working directory to the folder holding the
## raw WGS input data (expected folder layout: see README), e.g.
## setwd("/path/to/WGS_raw_data")

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(tibble)
})

strip_key <- function(x) {
  str_replace(basename(x), "\\.(vcf\\.gz|vcf|tsv\\.gz|tsv|txt)$", "")
}

snv_dir   <- "WGS_input/renamed/common/snv"
indel_dir <- "WGS_input/renamed/common/indel"
sv_dir    <- "WGS_input/renamed/common/sv"
cnv_dir   <- "WGS_input/renamed/common/cnv"
gtf_path  <- "WGS_input/gencode.v19.annotation.gtf"
meta_path <- "WGS_input/metadata/final_dna_metadata.csv"  # output of S2_dna_metadata.R

OUT_DIR <- "WGS_input/preprocessed"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

### Metadata
meta_raw <- fread(meta_path, data.table = FALSE)

sample_id   <- as.character(meta_raw$sample)
patient_id  <- as.character(meta_raw$patient)
tissue_type <- as.character(meta_raw$tissue_type)
gender      <- as.character(meta_raw$gender)

sub_col <- "tumor_type" 

meta <- tibble(
  sample_id   = sample_id,
  patient_id  = patient_id,
  subtype     = as.character(meta_raw[[sub_col]]),
  ploidy      = as.numeric(meta_raw$ploidy), 
  tcc         = as.numeric(meta_raw$tcc) ,
  tissue_type = tissue_type,
  gender      = gender
) 

desired_levels <- c("DDLPS","WD/DDLPS","WDLPS","PLPS","MLPS")
meta <- meta %>%
  mutate(subtype = factor(subtype, levels = desired_levels, ordered = TRUE)) %>%
  arrange(subtype, patient_id, sample_id)

fwrite(meta, file.path(OUT_DIR, "dna_metadata_ordered.tsv"), sep = "\t")

table(meta$subtype)

### Build manifest (4 types per sample)

snv_files   <- list.files(snv_dir,   pattern="\\.(vcf|tsv)(\\.gz)?$", full.names = TRUE)
indel_files <- list.files(indel_dir, pattern="\\.(vcf|tsv)(\\.gz)?$", full.names = TRUE)
sv_files    <- list.files(sv_dir,    pattern="\\.tsv(\\.gz)?$",       full.names = TRUE)
cnv_files   <- list.files(cnv_dir,   pattern="\\.(txt|tsv)(\\.gz)?$", full.names = TRUE)

snv_map   <- data.frame(sample_id=strip_key(snv_files),   snv_path=snv_files, stringsAsFactors=FALSE)
indel_map <- data.frame(sample_id=strip_key(indel_files), indel_path=indel_files, stringsAsFactors=FALSE)
sv_map    <- data.frame(sample_id=strip_key(sv_files),    sv_path=sv_files, stringsAsFactors=FALSE)
cnv_map   <- data.frame(sample_id=strip_key(cnv_files),   cnv_path=cnv_files, stringsAsFactors=FALSE)

manifest <- Reduce(function(x, y) full_join(x, y, by="sample_id"),
                   list(snv_map, indel_map, sv_map, cnv_map)) %>%
  filter(!is.na(snv_path) & !is.na(indel_path) & !is.na(sv_path) & !is.na(cnv_path))


### check the samples that can be used for subsequent analysis
sample_order <- meta %>%
  filter(sample_id %in% manifest$sample_id) %>%
 # arrange(subtype) %>%
  pull(sample_id) %>%
  unique()

samples_use <- c(sample_order, setdiff(manifest$sample_id, sample_order))
manifest <- manifest %>%
  mutate(sample_id = as.character(sample_id)) %>%
  arrange(factor(sample_id, levels = samples_use))


## Save manifest
fwrite(manifest, file.path(OUT_DIR, "WGS_data_manifest.tsv"), sep = "\t")

