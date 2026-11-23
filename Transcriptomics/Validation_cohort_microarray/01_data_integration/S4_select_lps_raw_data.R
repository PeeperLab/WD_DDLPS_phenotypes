
library(dplyr)

# set working directory
## ---- Working directory ----------------------------------------------------
## Before running, set the working directory to the folder holding the
## raw microarray input data (expected folder layout: see README), e.g.
## setwd("/path/to/microarray_raw_data")

gse_ids <- c("GSE21050", "GSE21122", "GSE30929", "GSE71118", "GSE159659")

# helper: copy & rename CEL files for one dataset
copy_cel_one_gse <- function(gse_id,
                             metadata_dir = "meta/lps_metadata",
                             raw_root     = "raw_data",
                             out_root     = "data/lps_data") {
  
  message("=== Processing ", gse_id, " ===")
  
  # 1) metadata
  meta_path <- file.path(metadata_dir, paste0(gse_id, "_lps_metadata.csv"))
  meta <- read.csv(meta_path, stringsAsFactors = FALSE)
  sample_ids <- meta$Series_sample_id %>% as.character() %>% unique() %>% na.omit()
  
  # 2) source folder (RAW extracted dir)
  source_folder <- file.path(raw_root, gse_id, paste0(gse_id, "_RAW"))
  # 3) target folder
  target_folder <- file.path(out_root, paste0(gse_id, "_LPS"))
  dir.create(target_folder, recursive = TRUE, showWarnings = FALSE)
  
  # 4) list candidate CEL/CEL.gz
  all_files <- list.files(
    source_folder,
    pattern = "\\.CEL(\\.gz)?$",
    full.names = FALSE,
    ignore.case = TRUE
  )
  
  copied_n <- 0L
  
  for (sid in sample_ids) {
    # compatible matching:
    # 1) sid.CEL / sid.CEL.gz
    # 2) sid_*.CEL / sid_*.CEL.gz
    pat_exact <- paste0("^", sid, "\\.CEL(\\.gz)?$")
    pat_pref  <- paste0("^", sid, "_.*\\.CEL(\\.gz)?$")
    
    hits <- all_files[grepl(pat_exact, all_files, ignore.case = TRUE)]
    if (length(hits) == 0) hits <- all_files[grepl(pat_pref, all_files, ignore.case = TRUE)]
    
    if (length(hits) == 0) {
      warning(sprintf("[%s] No CEL matched for %s", gse_id, sid))
      next
    }
    if (length(hits) > 1) {
      warning(sprintf("[%s] Multiple CEL matched for %s: %s ; take the first.",
                      gse_id, sid, paste(hits, collapse = ", ")))
    }
    
    source_file <- file.path(source_folder, hits[1])
    
    # keep original extension (.CEL or .CEL.gz), but rename prefix to sid only
    ext <- sub(paste0("^", sid, "(_.*)?"), "", hits[1], ignore.case = TRUE) # not reliable
    # easier: detect gz
    is_gz <- grepl("\\.gz$", hits[1], ignore.case = TRUE)
    target_name <- if (is_gz) paste0(sid, ".CEL.gz") else paste0(sid, ".CEL")
    target_file <- file.path(target_folder, target_name)
    
    ok <- file.copy(source_file, target_file, overwrite = TRUE)
    if (!ok) warning(sprintf("[%s] Failed to copy: %s -> %s", gse_id, source_file, target_file))
    if (ok) copied_n <- copied_n + 1L
  }
  
  message("Done: ", gse_id, " | Copied & renamed: ", copied_n,
          " | Output: ", target_folder)
  
  invisible(list(gse_id = gse_id, copied = copied_n, out = target_folder))
}

# run all datasets
results <- lapply(gse_ids, copy_cel_one_gse)
