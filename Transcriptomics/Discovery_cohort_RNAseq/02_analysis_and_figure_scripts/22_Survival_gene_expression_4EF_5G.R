# This script performs survival analyses by dividing the transcriptomic discovery cohort by the expression of a specific gene. 
# All adipocytic samples (cluster 6) are excluded from the analysis. 

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(survival)
library(survminer)
library(patchwork)

# set working directory
set.seed(1507)
## ---- Working directory ----------------------------------------------------
## Before running, set the working directory to the "Transcriptomics/Discovery_cohort_RNAseq" folder of
## this repository (the folder that holds Inputs/ and Outputs/), e.g.
## setwd("/path/to/WD_DDLPS_phenotypes/Transcriptomics/Discovery_cohort_RNAseq")

output_dir <- "Outputs/Figures/Survival_gene_expression_fig_4EF_5G"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Load data
survival_data_aspera <- read_csv("Inputs/aspera_rna_metadata_with_survival.csv")
survival_data_TCGA <- read_csv("Inputs/TCGA_sarc_dd_survival_data.csv")
metadata <- read_csv("Inputs/Metadata_discovery_WD_DD.csv")
exprs_data    <- readRDS("Inputs/Batch_corrected_data_discovery_WD_DD.rds")
clusters6 <- read.csv("Outputs/Files/Consensus_clustering/Consensus_clustering.csv", row.names = 1)

metadata <- metadata %>%
  left_join(clusters6, by = "sampleID")
colnames(metadata)[ncol(metadata)] <- "clusters6"

# Remove cluster 6 (Adipocytic samples)
metadata_filtered <- metadata %>% filter(clusters6 != 6)
exprs_data <- exprs_data[, colnames(exprs_data) %in% metadata_filtered$sampleID]

# Process TCGA survival data
survival_data_TCGA <- survival_data_TCGA %>%
  mutate(
    surv_time_days = coalesce(days_to_death, days_to_last_followup),
    status = ifelse(!is.na(days_to_death), 1, 0)
  )

# Fig. 4E (CD34), Fig. 4F (TGFBR3), Fig. 5G (PROM1)
for (gene_of_interest in c("CD34", "TGFBR3", "PROM1")) {

  if (!gene_of_interest %in% rownames(exprs_data)) {
    stop(paste("Gene", gene_of_interest, "not found in expression data"))
  }

  # Extract expression
  gene_expr_df <- data.frame(
    sampleID   = colnames(exprs_data),
    expression = as.numeric(exprs_data[gene_of_interest, ])
  )

  # Stratify: extreme quartiles (bottom 25% vs top 25%)
  quartiles <- quantile(gene_expr_df$expression, probs = c(0.25, 0.75), na.rm = TRUE)

  gene_expr_df <- gene_expr_df %>%
    mutate(expression_group = case_when(
      expression <= quartiles[1] ~ "Low (Q1)",
      expression >= quartiles[2] ~ "High (Q4)",
      TRUE ~ "Middle"
    )) %>%
    filter(expression_group != "Middle")

  # Merge with survival data
  merged_data_aspera <- gene_expr_df %>%
    inner_join(survival_data_aspera, by = c("sampleID" = "fastq_id"))

  merged_data_TCGA <- gene_expr_df %>%
    inner_join(survival_data_TCGA, by = "sampleID")

  common_cols <- intersect(names(merged_data_aspera), names(merged_data_TCGA))
  merged_data <- rbind(
    merged_data_aspera[, common_cols],
    merged_data_TCGA[, common_cols]
  )

  # Apply time cutoff
  time_cutoff <- 2500

  merged_data <- merged_data %>%
    mutate(
      surv_time_orig = surv_time_days,
      status         = ifelse(surv_time_orig > time_cutoff & status == 1, 0, status),
      surv_time_days = pmin(surv_time_orig, time_cutoff),
    
      # Convert to years
      surv_time_years = surv_time_days / 365.25
    )

  merged_data$expression_group <- factor(merged_data$expression_group, levels = c("Low (Q1)", "High (Q4)"))
  colors <- c("Low (Q1)" = "#80B1D2", "High (Q4)" = "darkred")

  # Survival analysis
  surv_obj <- Surv(time = merged_data$surv_time_years, event = merged_data$status)
  fit <- survfit(surv_obj ~ expression_group, data = merged_data, conf.type = "log-log")

  # Kaplan-Meier plot
  plot1 <- ggsurvplot(
    fit,
    data                  = merged_data,
    pval                  = TRUE,
    pval.method           = TRUE,
    pval.size             = 7,
    pval.method.size      = 7,
    conf.int              = FALSE,
    risk.table            = TRUE,
    risk.table.col        = "expression_group",
    risk.table.palette    = colors,
    legend.labs           = levels(merged_data$expression_group),
    legend.title          = paste(gene_of_interest, "expression"),
    xlab                  = "Time (years)",
    ylab                  = "Survival Probability",
    palette               = colors,
    legend                = "right",
    font.main             = c(18, "bold"),
    font.x                = c(18, "plain"),
    font.y                = c(18, "plain"),
    font.tickslab         = 18,
    font.legend           = c(18, "plain"),
    risk.table.fontsize   = 5.5,
    break.time.by         = 1,
    xlim                  = c(0, time_cutoff / 365.25),
    title                 = paste("Survival by", gene_of_interest, "expression")
  )

  combined <- plot1$plot / plot1$table +
    plot_layout(heights = c(3, 1))

  print(combined)

  # Sample size summary
  cat("\nSample sizes per group:\n")
  print(table(merged_data$expression_group))

  cat("\nEvents per group:\n")
  print(table(merged_data$expression_group, merged_data$status))

  ggsave(
    filename = file.path(output_dir, paste0(gene_of_interest, "_survival.pdf")),
    plot = combined,
    width = 8,
    height = 10,
    units = "in")
}
