set.seed(123)
## ---- Working directory ----------------------------------------------------
## Before running, set the working directory to the folder holding the
## raw WGS input data (expected folder layout: see README), e.g.
## setwd("/path/to/WGS_raw_data")

suppressPackageStartupMessages({
  library(stringr)
  library(dplyr)
  library(tibble)
})

### input dir
SNV_DIR   <- "WGS_input/snv/somatic_functional"
INDEL_DIR <- "WGS_input/indels/somatic_functional"
SV_DIR    <- "WGS_input/sophia/high_conf"
CNV_DIR   <- "WGS_input/aceseq_cnv"

### output dir for renamed files
OUT_DIR <- "WGS_input"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

OUT_SNV   <- file.path(OUT_DIR, "renamed", "snv")
OUT_INDEL <- file.path(OUT_DIR, "renamed", "indel")
OUT_SV    <- file.path(OUT_DIR, "renamed", "sv")
OUT_CNV   <- file.path(OUT_DIR, "renamed", "cnv_most_important_info")

dir.create(OUT_SNV,   showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_INDEL, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_SV,    showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_CNV,   showWarnings = FALSE, recursive = TRUE)


# Extract H021-xxxxx between 1st and 2nd "_"
extract_core_id <- function(fname) {
  m <- str_match(fname, "^[^_]+_(H021-[^_]+)_")
  m[, 2]
}

# Extract token after "__" (supports "__tumor03_blood" and "__tumor03__blood")
extract_double_underscore_token <- function(fname) {
  m1 <- str_match(fname, "__([^_]+?)_")  
  m1[,2]
}

# token -> suffix: tumor03 -> -T3; metastasis02 -> -M2; tumor -> -T1; metastasis -> -M1
token_to_suffix <- function(token) {
  if (is.na(token) || token == "") return("")
  token <- tolower(token)
  
  if (str_detect(token, "^tumor")) {
    num <- str_match(token, "^tumor(\\d+)$")[, 2]
    if (is.na(num)) num <- "1"
    num <- as.integer(num)
    return(paste0("-T", num))
  }
  
  if (str_detect(token, "^metastasis")) {
    num <- str_match(token, "^metastasis(\\d+)$")[, 2]
    if (is.na(num)) num <- "1"
    num <- as.integer(num)
    return(paste0("-M", num))
  }
  
  safe <- str_replace_all(token, "[^A-Za-z0-9]+", "")
  if (safe == "") safe <- "UNK"
  paste0("-X", safe)
}

# keep extension (supports .vcf.gz / .tsv.gz)
get_full_ext <- function(fname) {
  if (str_ends(fname, ".vcf.gz")) return("vcf.gz")
  if (str_ends(fname, ".tsv.gz")) return("tsv.gz")
  ext <- tools::file_ext(fname)
  if (ext == "") "" else ext
}

# Make destination names unique if collisions remain
make_unique_names <- function(x) {
  # x: character vector of file basenames
  # If duplicates exist, append -rep2, -rep3...
  base <- str_replace(x, "\\.(vcf\\.gz|tsv\\.gz|vcf|tsv|txt)$", "")
  ext  <- str_match(x, "\\.(vcf\\.gz|tsv\\.gz|vcf|tsv|txt)$")[,2]
  ext[is.na(ext)] <- ""
  
  out <- x
  dup_groups <- split(seq_along(x), base)
  for (idx in dup_groups) {
    if (length(idx) <= 1) next
    for (k in seq_along(idx)) {
      if (k == 1) next
      newb <- paste0(base[idx[k]], "-rep", k)
      out[idx[k]] <- if (ext[idx[k]] == "") newb else paste0(newb, ".", ext[idx[k]])
    }
  }
  out
}

# Generic rename+copy for SNV/INDEL/SV
rename_copy_folder <- function(in_dir, out_dir, pattern_regex, dry_run = TRUE) {
  files <- list.files(in_dir, full.names = TRUE, recursive = FALSE)
  files <- files[str_detect(basename(files), pattern_regex)]
  
  if (length(files) == 0) {
    message("[SKIP] no matched files in: ", in_dir)
    return(invisible(NULL))
  }
  
  df <- tibble(
    src_path = files,
    fname    = basename(files),
    core_id  = extract_core_id(basename(files)),
    token    = extract_double_underscore_token(basename(files)),
    ext      = vapply(basename(files), get_full_ext, character(1))
  )
  
  bad <- which(is.na(df$core_id))
  if (length(bad) > 0) {
    warning("Some files failed to extract core_id (H021-xxxxx):\n",
            paste(df$fname[bad], collapse = "\n"))
  }
  
  dup_ids <- df$core_id[duplicated(df$core_id) | duplicated(df$core_id, fromLast = TRUE)]
  
  df <- df %>%
    mutate(
      is_dup   = core_id %in% dup_ids,
      suffix   = ifelse(is_dup, vapply(token, token_to_suffix, character(1)), ""),
      new_base = paste0(core_id, suffix),
      new_name = ifelse(ext == "", new_base, paste0(new_base, ".", ext))
    )
  
  # ensure uniqueness if still collision
  df$new_name <- make_unique_names(df$new_name)
  df$dst_path <- file.path(out_dir, df$new_name)
  
  # preview
  print(df %>% select(fname, core_id, token, is_dup, new_name) %>% head(30))
  message("Total files: ", nrow(df), " | Duplicated IDs: ", length(unique(dup_ids)))
  
  if (!dry_run) {
    ok <- file.copy(df$src_path, df$dst_path, overwrite = TRUE)
    message("Copied: ", sum(ok), "/", length(ok), " to ", out_dir)
  } else {
    message("[DRY RUN] No files copied. Set dry_run=FALSE to execute.")
  }
  
  invisible(df)
}


# CNV: only _most_important_info
# Rule: keep only before "_most_important_info" + ".txt"

rename_copy_cnv_most_important <- function(in_dir, out_dir, dry_run = TRUE) {
  files <- list.files(in_dir, full.names = TRUE, recursive = FALSE)
  files <- files[str_detect(basename(files), "_most_important_info.*\\.txt$")]
  
  if (length(files) == 0) {
    message("[SKIP] no _most_important_info*.txt in: ", in_dir)
    return(invisible(NULL))
  }
  
  df <- tibble(
    src_path = files,
    fname    = basename(files),
    new_name = str_replace(basename(files), "_most_important_info.*\\.txt$", ".txt")
  )
  
  # ensure uniqueness
  df$new_name <- make_unique_names(df$new_name)
  df$dst_path <- file.path(out_dir, df$new_name)
  
  # collision check (should be rare)
  if (any(duplicated(df$dst_path))) {
    warning("Destination name collisions still detected in CNV after make_unique_names().")
  }
  
  print(df %>% select(fname, new_name) %>% head(30))
  message("Total CNV most_important_info files: ", nrow(df))
  
  if (!dry_run) {
    ok <- file.copy(df$src_path, df$dst_path, overwrite = TRUE)
    message("Copied: ", sum(ok), "/", length(ok), " to ", out_dir)
  } else {
    message("[DRY RUN] No files copied. Set dry_run=FALSE to execute.")
  }
  
  invisible(df)
}


# Run for 4 types

rename_copy_folder(INDEL_DIR, OUT_INDEL, pattern_regex = "\\.vcf(\\.gz)?$", dry_run = FALSE)
rename_copy_folder(SNV_DIR,   OUT_SNV,   pattern_regex = "\\.vcf(\\.gz)?$", dry_run = FALSE)
rename_copy_folder(SV_DIR,    OUT_SV,    pattern_regex = "\\.tsv(\\.gz)?$", dry_run = FALSE)

rename_copy_cnv_most_important(CNV_DIR, OUT_CNV, dry_run = FALSE)





### check how many ocerlapped samples among 4 types

get_sample_keys <- function(dir, pattern) {
  f <- list.files(dir, pattern = pattern, full.names = FALSE)
  ##remove .vcf / .vcf.gz / .tsv / .tsv.gz / .txt
  key <- str_replace(f, "\\.(vcf\\.gz|vcf|tsv\\.gz|tsv|txt)$", "")
  unique(key)
}

snv_keys   <- get_sample_keys(OUT_SNV,   "\\.(vcf|vcf\\.gz)$")
indel_keys <- get_sample_keys(OUT_INDEL, "\\.(vcf|vcf\\.gz)$")
sv_keys    <- get_sample_keys(OUT_SV,    "\\.(tsv|tsv\\.gz)$")
cnv_keys   <- get_sample_keys(OUT_CNV,   "\\.txt$")

cat("SNV samples  :", length(snv_keys), "\n")
cat("INDEL samples:", length(indel_keys), "\n")
cat("SV samples   :", length(sv_keys), "\n")
cat("CNV samples  :", length(cnv_keys), "\n")

common_4 <- Reduce(intersect, list(snv_keys, indel_keys, sv_keys, cnv_keys))
cat("\nCommon in all 4 folders:", length(common_4), "\n")


writeLines(common_4, "WGS_input/common_samples_all4.txt")





### copy common_4 files into one common folder (by type)

COMMON_DIR <- "WGS_input/renamed/common"

COMMON_SNV   <- file.path(COMMON_DIR, "snv")
COMMON_INDEL <- file.path(COMMON_DIR, "indel")
COMMON_SV    <- file.path(COMMON_DIR, "sv")
COMMON_CNV   <- file.path(COMMON_DIR, "cnv")

dir.create(COMMON_SNV,   showWarnings = FALSE, recursive = TRUE)
dir.create(COMMON_INDEL, showWarnings = FALSE, recursive = TRUE)
dir.create(COMMON_SV,    showWarnings = FALSE, recursive = TRUE)
dir.create(COMMON_CNV,   showWarnings = FALSE, recursive = TRUE)

# helper: copy files whose basename (w/o extension) is in common_4
copy_common_files <- function(src_dir, dst_dir, ext_regex, keys) {
  files <- list.files(src_dir, pattern = ext_regex, full.names = TRUE)
  if (length(files) == 0) {
    message("[SKIP] No files in ", src_dir)
    return(invisible(0))
  }
  
  base <- basename(files)
  file_key <- str_replace(base, "\\.(vcf\\.gz|vcf|tsv\\.gz|tsv|txt)$", "")
  keep <- file_key %in% keys
  
  files_keep <- files[keep]
  if (length(files_keep) == 0) {
    message("[SKIP] No common files matched in ", src_dir)
    return(invisible(0))
  }
  
  dst_paths <- file.path(dst_dir, basename(files_keep))
  ok <- file.copy(files_keep, dst_paths, overwrite = TRUE)
  message("Copied ", sum(ok), "/", length(ok), " -> ", dst_dir)
  invisible(sum(ok))
}

# copy each type
copy_common_files(OUT_SNV,   COMMON_SNV,   "\\.(vcf|vcf\\.gz)$", common_4)
copy_common_files(OUT_INDEL, COMMON_INDEL, "\\.(vcf|vcf\\.gz)$", common_4)
copy_common_files(OUT_SV,    COMMON_SV,    "\\.(tsv|tsv\\.gz)$", common_4)
copy_common_files(OUT_CNV,   COMMON_CNV,   "\\.txt$",            common_4)



