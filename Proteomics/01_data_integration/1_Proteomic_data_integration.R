# This script integrates different proteomic batches, performs quality controls, normalization and imputation,
# and outputs a combined dataset for downstream analysis.

library(readxl)
library(ggplot2)
library(data.table)
library(reshape2)
library(dplyr)
library(openxlsx)

set.seed(1507)
## Working directory 
## Before running, set the working directory to the "Proteomics" folder of
## this repository (the folder that holds Inputs/ and Outputs/), e.g.
## setwd("/path/to/WD_DDLPS_phenotypes/Proteomics")

#store downloaded data in Inputs folder and load:
batch968 <- read.xlsx("Inputs/Project_968_878_DIANN22_wo_cont_proteome_untransformed_export_proteins_wo_LPS18_TJ.xlsx")
batch1055 <- read.xlsx("Inputs/Project_1055_DIANN22_wo_cont_proteome_untransformed_export_proteins_TJ.xlsx")
batch1172 <- read.xlsx("Inputs/Project_1172_DIANN22_wo_cont_proteome_untransformed_export_proteins_TJ.xlsx")

# Ensure the output directory exists
output_dir <- "Inputs"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# batch968: melt and clean
expr_cols <- grep("^#", names(batch968), value = TRUE)

m_batch968 <- melt(setDT(batch968),
                   id.vars = names(batch968)[1],
                   measure.vars = expr_cols)

m_batch968 <- as.data.table(m_batch968)
m_batch968[, ID := sub(".*#\\d+_", "", variable)]
setnames(m_batch968, old = c("variable", "value"), new = c("Sample_Column", "value"))

# Basic cleanup
m_batch968 <- m_batch968[, .(Gene = get(names(batch968)[1]), ID, value)]
m_batch968[, value := as.numeric(value)]

# Log2 transform
m_batch968[, Log2_Value := log2(value + 1e-8)]
m_batch968[, Batch := "968"]

# 🔹 Collapse multiple isoforms per Gene per Sample by mean log2 value
m_batch968 <- m_batch968[, .(Log2_Value = mean(Log2_Value, na.rm = TRUE)),
                         by = .(Gene, ID, Batch)]


# batch1055: melt and clean
expr_cols <- names(batch1055)[-1]

m_batch1055 <- melt(setDT(batch1055),
                    id.vars = names(batch1055)[1],
                    measure.vars = names(batch1055)[-1])
setnames(m_batch1055, old = c(names(batch1055)[1], "variable", "value"),
         new = c("Gene", "ID", "value"))
m_batch1055 <- as.data.table(m_batch1055)
m_batch1055[, value := as.numeric(value)]
m_batch1055[, Log2_Value := log2(value + 1e-8)]
m_batch1055[, Batch := "1055"]
m_batch1055 <- m_batch1055[, .(Log2_Value = mean(Log2_Value, na.rm = TRUE)),
                           by = .(Gene, ID, Batch)]


# batch1172: melt and clean
expr_cols <- names(batch1172)[-1]

m_batch1172 <- melt(setDT(batch1172),
                    id.vars = names(batch1172)[1],
                    measure.vars = names(batch1172)[-1])
setnames(m_batch1172, old = c(names(batch1172)[1], "variable", "value"),
         new = c("Gene", "ID", "value"))
m_batch1172 <- as.data.table(m_batch1172)
m_batch1172[, value := as.numeric(value)]
m_batch1172[, Log2_Value := log2(value + 1e-8)]
m_batch1172[, Batch := "1172"]
m_batch1172 <- m_batch1172[, .(Log2_Value = mean(Log2_Value, na.rm = TRUE)),
                           by = .(Gene, ID, Batch)]


# Remove proteins with >75% missing values in each batch
filter_missing <- function(dt, value_col = "Log2_Value", threshold = 0.75) {
  # Calculate proportion of missing values per gene
  miss_stats <- dt[, .(missing_frac = mean(is.na(get(value_col)))), by = Gene]
  
  # Keep genes with <= threshold missing values
  keep_genes <- miss_stats[missing_frac <= threshold, Gene]
  
  message("Retaining ", length(keep_genes), " proteins (",
          round(100 * (1 - threshold), 0),
          "%+ detected) out of ", length(unique(dt$Gene)), " total.")
  
  dt[Gene %in% keep_genes]
}

# Apply to both batches
m_batch968 <- filter_missing(m_batch968, value_col = "Log2_Value", threshold = 0.75)
m_batch1055 <- filter_missing(m_batch1055, value_col = "Log2_Value", threshold = 0.75)
m_batch1172 <- filter_missing(m_batch1172, value_col = "Log2_Value", threshold = 0.75)

# Find common genes
common_genes <- intersect(
  intersect(unique(m_batch968$Gene), unique(m_batch1055$Gene)),
  unique(m_batch1172$Gene)
)

# Subset each batch to shared genes
m_batch968 <- m_batch968[Gene %in% common_genes]
m_batch1055 <- m_batch1055[Gene %in% common_genes]
m_batch1172 <- m_batch1172[Gene %in% common_genes]

# Merge batches
all_batches <- rbindlist(list(m_batch968, m_batch1055, m_batch1172), use.names = TRUE, fill = TRUE)

# Convert to linear
all_batches[, linear := 2^Log2_Value]

# Albumin abundance sample exclusion
plasma_proteins <- c("ALB")

# Compute totals and global median
totals_all <- all_batches[, .(total_linear = sum(linear, na.rm = TRUE)), by = ID]
global_total_all <- median(totals_all$total_linear, na.rm = TRUE)
all_batches <- merge(all_batches, totals_all, by = "ID")

# Compute total plasma abundance per sample
plasma_abundance <- all_batches[Gene %in% plasma_proteins,
                                .(plasma_total = sum(linear, na.rm = TRUE)),
                                by = ID]

# Join plasma totals with sample totals
plasma_abundance <- merge(plasma_abundance,
                          unique(all_batches[, .(ID, total_linear)]),
                          by = "ID")

# Calculate fraction
plasma_abundance[, plasma_fraction := plasma_total / total_linear]

# Albumin abundance sample exclusion
bad_samples <- plasma_abundance[plasma_fraction > 0.25, ID]

all_batches_filtered <- all_batches[!ID %in% bad_samples]
all_batches <- all_batches_filtered

# Align all batches using a common global total reference
# Calculate scaling factor and apply
sample_totals <- all_batches[, .(total_linear = sum(linear, na.rm = TRUE)), by = ID]
global_median <- median(sample_totals$total_linear, na.rm = TRUE)

all_batches <- merge(all_batches, sample_totals, by = "ID", suffixes = c("", "_new"))
all_batches[, linear := linear / total_linear_new * global_median]
all_batches[, total_linear_new := NULL]

# Convert back to log2 for imputation
all_batches[, Log2_Value := log2(linear + 1e-8)]
all_batches[, c("linear", "total_linear") := NULL]

# Split batches
m_batch968  <- all_batches[Batch == "968"]
m_batch1055 <- all_batches[Batch == "1055"]
m_batch1172 <- all_batches[Batch == "1172"]

# Batch-wise imputation
impute_perseus <- function(dt, log2_col = "Log2_Value", batch_name = NULL, width = 0.3, downshift = 2.4, seed = 1507){
  set.seed(seed)
  dt[, (log2_col) := as.numeric(get(log2_col))]
  imputed_col <- paste0(log2_col, "_Imputed")
  dt[, (imputed_col) := get(log2_col)]
  dt[, (imputed_col) := {
    x <- get(imputed_col)
    na_idx <- which(is.na(x))
    if(length(na_idx) > 0){
      mu <- mean(x, na.rm = TRUE)
      sigma <- sd(x, na.rm = TRUE)
      x[na_idx] <- rnorm(length(na_idx), mean = mu - downshift*sigma, sd = width*sigma)
    }
    x
  }, by = ID]
  if(!is.null(batch_name)) dt[, Batch := batch_name]
  return(dt)
}

m_batch968 <- impute_perseus(m_batch968, batch_name = "968")
m_batch1055 <- impute_perseus(m_batch1055, batch_name = "1055")
m_batch1172 <- impute_perseus(m_batch1172, batch_name = "1172")

combined <- rbindlist(
  list(
    m_batch968[, .(Gene, ID, Log2_Value, Log2_Value_Imputed, Batch)],
    m_batch1055[, .(Gene, ID, Log2_Value, Log2_Value_Imputed, Batch)],
    m_batch1172[, .(Gene, ID, Log2_Value, Log2_Value_Imputed, Batch)]
  ),
  use.names = TRUE,
  fill = TRUE
)
combined <- as.data.table(combined)

# Convert to linear
combined[, linear_Imputed := 2^Log2_Value_Imputed]

# Pivot to wide format: Genes as rows, samples as columns
wide <- dcast(combined, ID ~ Gene, value.var = "linear_Imputed")

# Write to Excel
write.xlsx(wide, file = file.path(output_dir, "Proteomics_Integrated.xlsx"), rowNames = FALSE)

# Filter metadata to contain the same samples (if samples are excluded)
metadata <- read.xlsx("Inputs/Proteomic_metadata.xlsx")
filtered_metadata <- metadata[metadata$ID %in% wide$ID, ]

# Write to Excel
write.xlsx(filtered_metadata, file = file.path(output_dir, "Proteomics_metadata.xlsx"), rowNames = FALSE)

