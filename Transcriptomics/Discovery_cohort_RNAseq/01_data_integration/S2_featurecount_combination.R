
### Batch-preprocess featureCounts for ALL datasets under raw_featurecount/
### - Auto-detect datasets (subfolders)
### - Auto-match metadata: metadata/filtered_SraRunTable_<DATASET>.txt
### - Save:
###     rds/counts_data_<TAG>.rds        (TAG = 282597 for PRJNA282597; "german" stays "german")
###     metadata/metadata_<TAG>.csv      (skip german)
### - Convert Ensembl -> Gene Symbol (biomaRt), drop rowSum=0, sum duplicate symbols

###Load the necessary libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(biomaRt)
})

### Set working directory
## ---- Working directory ----------------------------------------------------
## Before running, set the working directory to the folder holding the
## raw RNA-seq input data (expected folder layout: see README), e.g.
## setwd("/path/to/RNAseq_raw_data")

### Define root directories and create output folders
counts_root <- "data/raw_featurecount"
meta_root   <- "meta/metadata" 
### out dir
out_dir1   <- "meta/metadata" 
out_dir2 <- "data/rds"
for (d in c(out_dir1, out_dir2)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

### Connect to Ensembl 
ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

### Discover all datasets 
datasets <- list.dirs(counts_root, full.names = FALSE, recursive = FALSE)
datasets <- datasets[datasets != "" & !startsWith(datasets, ".")]


### Loop through each dataset
for (ds in datasets) {
  
  folder_path <- file.path(counts_root, ds)
  
  file_list <- list.files(path = folder_path, pattern = "\\.biotype_counts_mqc\\.(tsv|txt)$", full.names = TRUE)
  
  ### Read and Merge Count Files
  merged_data <- NULL
  
  for (file in file_list) {
    # Extract sample name
    sample_name <- gsub("\\.biotype_counts_mqc\\.(tsv|txt)$", "", basename(file))
    
    # Read file and select rows starting with ENSG
    data <- read.table(file, comment.char = "#", header = FALSE, stringsAsFactors = FALSE)
    data <- data[grep("^ENSG", data$V1), ]
    colnames(data) <- c("ensembl_gene_id", sample_name)
    
    # Merge
    if (is.null(merged_data)) {
      merged_data <- data
    } else {
      merged_data <- merge(merged_data, data, by = "ensembl_gene_id", all = TRUE)
    }
  }
  
  # Replace any NA values generated during merge with 0
  merged_data[is.na(merged_data)] <- 0
  
  ###: Convert Ensembl ID to Gene Symbol
  ensembl_gene_ids <- merged_data$ensembl_gene_id
  
  gene_info <- getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol"),
    filters    = "ensembl_gene_id",
    values     = ensembl_gene_ids,
    mart       = ensembl
  )
  
  colnames(gene_info)[colnames(gene_info) == "hgnc_symbol"] <- "gene_symbol"
  raw_count <- inner_join(gene_info, merged_data, by = "ensembl_gene_id")
  
  ### Pre-process Count Data
  # 1. Remove rows without gene symbol
  raw_count_noNAs <- raw_count[raw_count$gene_symbol != "", ]
  
  # 2. Delete rows where sum of all expression values is 0
  # (Columns 3 to end are the numeric sample columns)
  row_sums <- rowSums(raw_count_noNAs[, 3:ncol(raw_count_noNAs)])
  nonsum0_count_data <- raw_count_noNAs[row_sums != 0, ]
  
  # 3. Sum duplicate genes
  # using where(is.numeric) is safer than starts_with("SRR") in case you have ERR/DRR samples
  count_data_duplicate_summed <- nonsum0_count_data %>%
    group_by(gene_symbol) %>%
    summarise(across(where(is.numeric), sum, na.rm = TRUE))
  
  # Save RDS
  ds_tag <- sub("^PRJNA", "", ds) # Format to purely numeric tag if it starts with PRJNA
  out_rds <- file.path(out_dir2, paste0("counts_data_", ds_tag, ".rds"))
  saveRDS(count_data_duplicate_summed, file = out_rds)
  
  ### Handle Metadata
  if (tolower(ds) == "german") {
    message("-> Skipping metadata for 'german' dataset as requested.")
    next
  }
  
  meta_candidates <- c(
    file.path(meta_root, paste0("filtered_SraRunTable_", ds, ".csv")),
    file.path(meta_root, paste0("filtered_SraRunTable_PRJNA", ds_tag, ".csv"))
  )
  meta_path <- meta_candidates[file.exists(meta_candidates)][1]
  
  raw_metadata <- read_csv(meta_path, show_col_types = FALSE)
  
  #  Customizable Subtype Extraction based on Dataset 
  if (grepl("282597", ds)) {
    selected_metadata <- data.frame(sampleID = raw_metadata$Run, subtype = raw_metadata$disease)
    selected_metadata$subtype[selected_metadata$subtype == "Liposarcoma - dedifferentiated"] <- "DDLPS"
    selected_metadata$subtype[selected_metadata$subtype == "Liposarcoma - pleomorphic"] <- "pleomorphic"
    
  } else if (grepl("987378", ds)) {
    selected_metadata <- data.frame(sampleID = raw_metadata$Run, subtype = raw_metadata$Histologic_subtypes)
    
  } else if (grepl("1053633", ds)) {
    selected_metadata <- data.frame(sampleID = raw_metadata$Run, subtype = raw_metadata$tissue)
    selected_metadata$subtype <- substr(selected_metadata$subtype, 1, 5)
    
  } else {
    message("   No specific metadata mapping defined for ", ds, ". Using 'Run' and entirely keeping original.")
    if ("Run" %in% colnames(raw_metadata)) {
      selected_metadata <- data.frame(sampleID = raw_metadata$Run, subtype = "Unknown")
    } else {
      selected_metadata <- data.frame(sampleID = raw_metadata[[1]], subtype = "Unknown")
    }
  }
  
  # Add dataset name and reorder to match counts data columns
  selected_metadata$dataset <- ds
  
  # Only keep samples that actually exist in the count data, and order them
  count_samples <- colnames(count_data_duplicate_summed)[-1] # skip 'gene_symbol'
  selected_metadata <- selected_metadata %>%
    filter(sampleID %in% count_samples) %>%
    arrange(match(sampleID, count_samples))
  
  # Save Metadata
  out_csv <- file.path(out_dir1, paste0("metadata_", ds_tag, ".csv"))
  write_csv(selected_metadata, file = out_csv)
}

message("\nAll datasets processed successfully!")