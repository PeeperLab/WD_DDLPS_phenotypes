
### This code is used to deal with the metadata of fastq datasets. 
###It can be used to select the liposarcoma and normal fat samples from 
###the BioProject's metadata and can be used to generate the srr_id list.

# set working directory
## ---- Working directory ----------------------------------------------------
## Before running, set the working directory to the folder holding the
## raw RNA-seq input data (expected folder layout: see README), e.g.
## setwd("/path/to/RNAseq_raw_data")
# input and output file path
input_path  <- "meta/raw_metadata/NCBI"   
output_dir <- "meta/raw_metadata/NCBI/output" 
output_dir2 <- "meta/metadata" 
for (d in c(output_dir, output_dir2)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# get all the metadata.csv files
csv_files <- list.files(input_path, pattern = "\\.csv$")

for (file in csv_files) {
  
  full_file_path <- file.path(input_path, file)
  data <- read.csv(full_file_path, stringsAsFactors = FALSE)
  colnames(data) <- gsub(" ", ".", colnames(data))
  
  filtered_data <- data.frame()
  if (basename(file) == "SraRunTable_PRJNA1053633.csv") {
    filtered_data <- data
  } else if ("disease" %in% colnames(data)) {
    filtered_data <- subset(data, Assay.Type == "RNA-Seq" & grepl("liposarcoma|Liposar", disease, ignore.case = TRUE))
  } else if ("isolate" %in% colnames(data)) {
    filtered_data <- subset(data, Assay.Type == "RNA-Seq" & grepl("liposarcoma|Liposar", isolate, ignore.case = TRUE))
  } else if ("tissue" %in% colnames(data)) {
    filtered_data <- subset(data, grepl("RNA", Assay.Type, ignore.case = TRUE) &
                              grepl("(DDLPS\\s*\\d+\\s*tumou?r\\s*tissue|WDLPS\\s*\\d+\\s*tumou?r\\s*tissue|liposarcoma|\\bLPS\\b|peritumou?ral[ _-]?fat|adipose|normal[ _-]?fat)",
                                    tissue, ignore.case = TRUE, perl = TRUE))
  } 
  # write result if data found
  if (nrow(filtered_data) == 0) {
    print(paste("No matching records found in", file, "based on the given conditions."))
  } else {
    # output filtered csv
    filtered_csv_file <- paste0("filtered_", file)
    full_filtered_csv_file_path <- file.path(output_dir2, filtered_csv_file)
    write.csv(filtered_data, full_filtered_csv_file_path, row.names = FALSE)
    
    # output SRR IDs for fastq downloading
    srr_id <- filtered_data$Run
    acclist_txt_file <- sub("SraRunTable", "acclist", file)
    acclist_txt_file <- sub("\\.csv$", "_srr_ids.txt", acclist_txt_file)
    full_acclist_txt_file_path <- file.path(output_dir, acclist_txt_file)
    write.table(srr_id, full_acclist_txt_file_path, row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
}
