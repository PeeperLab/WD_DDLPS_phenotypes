### This script is used to combine the metadata of all microarray datasets.

library(dplyr)
library(readr)

# set working directory
## ---- Working directory ----------------------------------------------------
## Before running, set the working directory to the folder holding the
## raw microarray input data (expected folder layout: see README), e.g.
## setwd("/path/to/microarray_raw_data")

input_dir <- "meta/lps_metadata"
out_dir <- "meta/metadata"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# list all *_lps_metadata.csv files
csv_files <- list.files(
  path = input_dir,
  pattern = "_lps_metadata\\.csv$",
  full.names = TRUE
)

# dataset -> microarray platform mapping
platform_map <- c(
  GSE30929  = "HG-U133A",
  GSE21122  = "HG-U133A",
  GSE21050  = "HG-U133_Plus_2",
  GSE71118  = "HG-U133_Plus_2",
  GSE159659 = "Clariom_S_Human"
)

all_list <- list()

for (f in csv_files) {
  # dataset name = remove "_lps_metadata.csv"
  dataset <- sub("_lps_metadata\\.csv$", "", basename(f))
  
  data <- read_csv(f, show_col_types = FALSE)
  selected <- data %>%
    dplyr::select(Series_sample_id, subtype) %>%
    mutate(
      dataset = dataset,
      microarray = unname(platform_map[dataset])
    )
  
  # dataset-specific fix
  if (dataset == "GSE71118") {
    selected <- selected %>%
      mutate(subtype = ifelse(subtype == "Liposarcoma", "Dedifferentiated liposarcoma", subtype))
  }
  
  all_list[[dataset]] <- selected
}

integrate_metadata <- bind_rows(all_list)

# normalize subtype strings
integrate_metadata <- integrate_metadata %>%
  mutate(subtype = tolower(as.character(subtype)))

# change subtype names to consistent names
integrate_metadata <- integrate_metadata %>%
  mutate(
    subtype = case_when(
      subtype %in% c("myxoid/round cell", "myxoid/rc", "myxoid") ~ "MLPS",
      subtype %in% c("dedifferentiated liposarcoma", "dedifferentiated") ~ "DDLPS",
      subtype %in% c("well differentiated liposarcoma", "well-differentiated", "well differentiated") ~ "WDLPS",
      subtype %in% c("adipose tissue", "normal_control", "normal control", NA) ~ "Fat_tissue",
      subtype%in% c("pleomorphic") ~ "PLPS",
      TRUE ~ subtype
    )
  )

## exclude the low quality samples
low_qc_sample <- c("GSM528276", "GSM766533", "GSM1827111", "GSM4837054", "GSM4837069")

selected_data <- integrate_metadata[!integrate_metadata$Series_sample_id %in% low_qc_sample, ]

colnames(selected_data)[1] <- "sampleID" 

# Dataset labels used in the manuscript (BioProject accession numbers)
dataset_labels <- c(GSE30929  = "PRJNA146279",
                    GSE21122  = "PRJNA129707",
                    GSE21050  = "PRJNA126435",
                    GSE71118  = "PRJNA290394",
                    GSE159659 = "PRJNA670101")
selected_data$dataset <- unname(dataset_labels[selected_data$dataset])
stopifnot(!anyNA(selected_data$dataset))
table(selected_data$dataset)

write_csv(selected_data, file.path(out_dir,"Metadata_validation_all.csv"))

