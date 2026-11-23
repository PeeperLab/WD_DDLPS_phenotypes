
### This script extracts liposarcoma samples from metadata CSV files.
library(readr)

## ---- Working directory ----------------------------------------------------
## Before running, set the working directory to the folder holding the
## raw microarray input data (expected folder layout: see README), e.g.
## setwd("/path/to/microarray_raw_data")

input_dir <- "meta/raw_metadata"
out_dir   <- "meta/lps_metadata"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# list all *_metadata.csv files
csv_files <- list.files(
  path = input_dir,
  pattern = "_metadata\\.csv$",
  full.names = TRUE
)

for (f in csv_files) {
  # suffix: remove "_metadata.csv" from filename
  suffix <- sub("_metadata\\.csv$", "", basename(f))
  
  # read
  data <- read_csv(f, show_col_types = FALSE)
  
  # filter liposarcoma + adipose tissue + normal control
  keep_patterns <- c("liposarcoma", "adipose tissue", "normal control")
  
  liposarcoma_data <- subset(
    data,
    grepl(paste(keep_patterns, collapse = "|"), tolower(tumour_type))
  )
  
  # write
  out_file <- file.path(out_dir, paste0(suffix, "_lps_metadata.csv"))
  write.csv(liposarcoma_data, out_file, row.names = FALSE)
  
  message("Saved: ", out_file, " (n = ", nrow(liposarcoma_data), ")")
}

