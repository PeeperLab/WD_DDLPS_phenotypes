# This script combines all proteomic batches without excluding samples based on quality controls. No imputation is performed
# in order to analyse whether a protein is detected or not in a sample. 

library(data.table)
library(openxlsx)

set.seed(1507)
## ---- Working directory ----------------------------------------------------
## Before running, set the working directory to the "Proteomics" folder of
## this repository (the folder that holds Inputs/ and Outputs/), e.g.
## setwd("/path/to/WD_DDLPS_phenotypes/Proteomics")

output_dir <- "Inputs"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Load batches
batch968 <- read.xlsx("./Inputs/Project_968_878_DIANN22_wo_cont_proteome_untransformed_export_proteins_wo_LPS18_TJ.xlsx")
batch1055 <- read.xlsx("./Inputs/Project_1055_DIANN22_wo_cont_proteome_untransformed_export_proteins_TJ.xlsx")
batch1172 <- read.xlsx("./Inputs/Project_1172_DIANN22_wo_cont_proteome_untransformed_export_proteins_TJ.xlsx")

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

# Merge batches
combined <- rbindlist(
  list(m_batch968, m_batch1055, m_batch1172),
  use.names = TRUE,
  fill = TRUE
)

combined_wide <- dcast(combined, ID + Batch ~ Gene, value.var = "Log2_Value")
write.xlsx(combined_wide, file = file.path(output_dir, "Raw_combined_data.xlsx"), rowNames = FALSE)
