###This script is used to deal with the microarray metadata

# set working directory
## ---- Working directory ----------------------------------------------------
## Before running, set the working directory to the folder holding the
## raw microarray input data (expected folder layout: see README), e.g.
## setwd("/path/to/microarray_raw_data")

file_path <- "meta/raw_metadata"

# Function to extract sample id and other fields
extract_sample_ids <- function(content) {
  sample_id_line <- grep("!Series_sample_id", content, value = TRUE)
  sample_ids <- strsplit(sample_id_line, '"')[[1]][2]
  sample_ids <- unlist(strsplit(sample_ids, " "))
  return(sample_ids)
}

extract_field_multiline <- function(field_name, content) {
  start_line <- grep(field_name, content)
  
  field_lines <- content[start_line]
  for (i in (start_line + 1):length(content)) {
    if (grepl("^!", content[i])) break
    field_lines <- paste(field_lines, content[i], sep = " ")
  }
  
  field_values <- strsplit(field_lines, '"')[[1]]
  field_values <- field_values[seq(2, length(field_values), 2)]
  return(field_values)
}

# GSE30929, GSE21122, GSE21050: tumour type / subtype parsed from the sample characteristics
# (output includes a 'tissue' column)
process_file_characteristics <- function(file_path, output_path) {
  # read all the contents
  file_content <- readLines(file_path)
  
  # extract sample id
  sample_ids <- extract_sample_ids(file_content)
  
  sample_titles <- extract_field_multiline("!Sample_title", file_content)
  sample_status <- extract_field_multiline("!Sample_status", file_content)
  sample_submission_dates <- extract_field_multiline("!Sample_submission_date", file_content)
  sample_last_update_dates <- extract_field_multiline("!Sample_last_update_date", file_content)
  sample_types <- extract_field_multiline("!Sample_type", file_content)
  sample_platform_ids <- extract_field_multiline("!Sample_platform_id", file_content)
  sample_contact_names <- extract_field_multiline("!Sample_contact_name", file_content)
  sample_contact_institutes <- extract_field_multiline("!Sample_contact_institute", file_content)
  sample_extract_protocol_ch1s <- extract_field_multiline("!Sample_extract_protocol_ch1", file_content)
  
  tissue_info <- extract_field_multiline("!Sample_characteristics_ch1\t\"tissue:", file_content)
  tissue_values <- sapply(tissue_info, function(x) sub(".*tissue: ", "", x))
  
  if (grepl("GSE21122", file_path)) {
    disease_status_info <- extract_field_multiline("!Sample_characteristics_ch1\t\"disease status:", file_content)
    disease_status_values <- sapply(disease_status_info, function(x) {
      parts <- unlist(strsplit(sub(".*disease status: ", "", x), ":"))
      c(tumour_type = parts[1], subtype = ifelse(length(parts) > 1, parts[2], NA))
    })
    tumour_type <- disease_status_values["tumour_type", ]
    subtype <- disease_status_values["subtype", ]
  } else if (grepl("GSE30929", file_path)) {
    sample_source_names <- extract_field_multiline("!Sample_source_name_ch1", file_content)
    tumour_type <- sub(",.*", "", sample_source_names)
    subtype <- sub(".*, ", "", sample_source_names)
  } else if (grepl("GSE21050", file_path)) {
    disease_diagnosis_info <- extract_field_multiline("!Sample_characteristics_ch1\t\"diagnosis:", file_content)
    tumour_type <- vector("character", length(disease_diagnosis_info))
    subtype <- vector("character", length(disease_diagnosis_info))
    for (i in seq_along(disease_diagnosis_info)) {
      diagnosis <- sub(".*diagnosis: ", "", disease_diagnosis_info[i])
      if (grepl(" - ", diagnosis)) {
        parts <- unlist(strsplit(diagnosis, " - "))
        tumour_type[i] <- parts[1]
        subtype[i] <- parts[2]
      } else {
        tumour_type[i] <- diagnosis
        subtype[i] <- NA
      }
    }
  }
  
  df <- data.frame(
    Series_sample_id = sample_ids,
    Sample_title = sample_titles,
    Sample_status = sample_status,
    Sample_submission_date = sample_submission_dates,
    Sample_last_update_date = sample_last_update_dates,
    tissue = tissue_values,
    tumour_type = tumour_type,
    subtype = subtype,
    Sample_type = sample_types,
    Sample_platform_id = sample_platform_ids,
    Sample_contact_name = sample_contact_names,
    Sample_contact_institute = sample_contact_institutes,
    sample_extract_protocol_ch1 = sample_extract_protocol_ch1s,
    stringsAsFactors = FALSE
  )
  
  write.csv(df, output_path, row.names = FALSE, quote = TRUE)
}

# GSE71118, GSE159659: tumour type / subtype taken from the sample source name
# (no 'tissue' field; all fields recycled to the number of samples)
process_file_source_name <- function(file_path, output_path) {
  # read all the contents
  file_content <- readLines(file_path)
  
  # extract sample id
  sample_ids <- extract_sample_ids(file_content)
  
  sample_titles <- extract_field_multiline("!Sample_title", file_content)
  sample_status <- extract_field_multiline("!Sample_status", file_content)
  sample_submission_dates <- extract_field_multiline("!Sample_submission_date", file_content)
  sample_last_update_dates <- extract_field_multiline("!Sample_last_update_date", file_content)
  sample_types <- extract_field_multiline("!Sample_type", file_content)
  sample_platform_ids <- extract_field_multiline("!Sample_platform_id", file_content)
  sample_contact_names <- extract_field_multiline("!Sample_contact_name", file_content)
  sample_contact_institutes <- extract_field_multiline("!Sample_contact_institute", file_content)
  sample_extract_protocol_ch1s <- extract_field_multiline("!Sample_extract_protocol_ch1", file_content)
  tumour_type <- extract_field_multiline("!Sample_source_name_ch1", file_content)
  subtype <- extract_field_multiline("!Sample_source_name_ch1", file_content)
  
  
  # Ensure all extracted fields have the same length
  len <- length(sample_ids)
  sample_titles <- rep_len(sample_titles, len)
  sample_status <- rep_len(sample_status, len)
  sample_submission_dates <- rep_len(sample_submission_dates, len)
  sample_last_update_dates <- rep_len(sample_last_update_dates, len)
  tumour_type <- rep_len(tumour_type, len)
  subtype <- rep_len(subtype, len)
  sample_types <- rep_len(sample_types, len)
  sample_platform_ids <- rep_len(sample_platform_ids, len)
  sample_contact_names <- rep_len(sample_contact_names, len)
  sample_contact_institutes <- rep_len(sample_contact_institutes, len)
  sample_extract_protocol_ch1s <- rep_len(sample_extract_protocol_ch1s, len)
  
  # Create a data frame
  df <- data.frame(
    Series_sample_id = sample_ids,
    Sample_title = sample_titles,
    Sample_status = sample_status,
    Sample_submission_date = sample_submission_dates,
    Sample_last_update_date = sample_last_update_dates,
    tumour_type = tumour_type,
    subtype = subtype,
    Sample_type = sample_types,
    Sample_platform_id = sample_platform_ids,
    Sample_contact_name = sample_contact_names,
    Sample_contact_institute = sample_contact_institutes,
    sample_extract_protocol_ch1 = sample_extract_protocol_ch1s,
    stringsAsFactors = FALSE
  )
  
  # Write the data frame to the output file
  write.csv(df, output_path, row.names = FALSE, quote = TRUE)
}

# Datasets and the parser used for each
datasets <- c(
  GSE30929  = "characteristics",
  GSE21122  = "characteristics",
  GSE21050  = "characteristics",
  GSE71118  = "source_name",
  GSE159659 = "source_name"
)

# Process each file
for (gse in names(datasets)) {
  in_file  <- file.path(file_path, paste0(gse, "_series_matrix.txt"))
  out_file <- file.path(file_path, paste0(gse, "_metadata.csv"))
  if (datasets[[gse]] == "characteristics") {
    process_file_characteristics(in_file, out_file)
  } else {
    process_file_source_name(in_file, out_file)
  }
}
