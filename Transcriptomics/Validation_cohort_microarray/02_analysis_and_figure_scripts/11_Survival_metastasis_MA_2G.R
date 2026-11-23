# This script analyzes survival data for the 6 clusters identified in the consensus clustering analysis. 

library(readr)
library(tidyverse)
library(ggplot2)
library(survival)
library(survminer)
library(ggpubr)
library(patchwork)

# set working directory
set.seed(1507)
## ---- Working directory ----------------------------------------------------
## Before running, set the working directory to the "Transcriptomics/Validation_cohort_microarray" folder of
## this repository (the folder that holds Inputs/ and Outputs/), e.g.
## setwd("/path/to/WD_DDLPS_phenotypes/Transcriptomics/Validation_cohort_microarray")

output_dir <- "Outputs/Figures/Survival_clusters_2G"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

survival_data1 <- read_csv("Inputs/survival_data_GSE21050.csv")
survival_data2 <- read_csv("Inputs/survival_data_GSE71118.csv")
cluster_data <- read.csv("Outputs/Files/Consensus_clustering/Consensus_clustering.csv", row.names = 1)

survival_data <- rbind(survival_data1,survival_data2)
colnames(survival_data)[colnames(survival_data) == "sampleid"] <- "sampleID"
colnames(survival_data)[colnames(survival_data) == "time"] <- "surv_time_years"
survival_data$status <- ifelse(survival_data$metastasis == "yes", 1, 0)

cluster_data$consensus_clusters <- as.integer(cluster_data$consensus_clusters)

merged_df <- survival_data %>%
  left_join(cluster_data, by = "sampleID")

merged_df$consensus_clusters <- factor(
  paste0("Cluster_", merged_df$consensus_clusters),
  levels = paste0("Cluster_", 1:6)
)

colors <- c(
  "Cluster_1" = "#A5C1A4",
  "Cluster_2" = "#827BB3",
  "Cluster_3" = "#DF9677",
  "Cluster_4" = "#1F70B5",
  "Cluster_5" = "#F7CBDF",
  "Cluster_6" = "#FDEDB9"
)

# Apply time cutoff (e.g., 10 years)
time_cutoff <- 5  # Change this value as needed

merged_df <- merged_df %>%
  mutate(
    surv_time_orig = surv_time_years,
    
    # Censor events occurring after cutoff
    status = ifelse(surv_time_orig > time_cutoff & status == 1, 0, status),
    
    # Truncate time
    surv_time_years = pmin(surv_time_orig, time_cutoff)
  )

# Surv object
surv_obj <- Surv(time = merged_df$surv_time_years, event = merged_df$status)

fit <- survfit(surv_obj ~ consensus_clusters, data = merged_df, conf.type = "log-log")

# Kaplan-Meier curive
plot1 <- ggsurvplot(
  fit,
  data = merged_df,
  pval = TRUE,
  pval.method = TRUE,
  pval.size = 7,
  pval.method.size = 7,
  conf.int = FALSE,
  risk.table = TRUE,
  risk.table.col = "consensus_clusters",
  risk.table.palette = colors,
  legend.labs = levels(merged_df$consensus_clusters),
  legend.title = "Cluster",
  xlab = "Time (years)",
  ylab = "Survival Probability",
  palette = colors,
  legend = "right",
  font.main = c(18, "bold"),
  font.x = c(18, "plain"),
  font.y = c(18, "plain"),
  font.tickslab = 18,
  font.legend = c(18, "plain"),
  font.subtitle = c(18, "plain"),
  font.caption = c(18, "plain"),
  risk.table.fontsize = 5.5,
  xlim = c(0, time_cutoff),
)

print(plot1)

combined <- plot1$plot / plot1$table +
  plot_layout(heights = c(3, 1))

ggsave(
  filename = file.path(output_dir, "Survival_clusters.pdf"),
  plot = combined,
  width = 7.5,
  height = 6,
  units = "in")


pairwise_lr <- pairwise_survdiff(
  Surv(surv_time_years, status) ~ consensus_clusters,
  data = merged_df,
  p.adjust.method = "BH"
)

# Extract lower-triangle p-values
pval_lower <- pairwise_lr$p.value

clusters <- levels(merged_df$consensus_clusters)

# Initialize empty matrix
pval_matrix <- matrix(
  NA,
  nrow = length(clusters),
  ncol = length(clusters),
  dimnames = list(clusters, clusters)
)

# Fill lower triangle
pval_matrix[lower.tri(pval_matrix)] <- as.vector(pval_lower)

# Mirror to upper triangle
pval_matrix <- t(pval_matrix)
pval_matrix[lower.tri(pval_matrix)] <- as.vector(pval_lower)

# Diagonal = 1 (or NA if you prefer)
diag(pval_matrix) <- 1

pval_matrix_round <- signif(pval_matrix, 3)
pval_matrix_round

# Extract survival table from survfit object
sf_df <- surv_summary(fit, data = merged_df)

# Find time point where n.risk <= 5 per cluster
cutoff_per_cluster <- sf_df %>%
  group_by(strata) %>%
  filter(n.risk <= 5) %>%
  summarise(
    cutoff_time = min(time),
    .groups = "drop"
  )

# View results
print(cutoff_per_cluster)

# Cluster 2 vs all other clusters combined (excl. cluster 4)
merged_data_binary <- merged_df %>%
  #filter(consensus_clusters != "Cluster_4") %>%
  mutate(cluster_group = ifelse(consensus_clusters == "Cluster_2", "Cluster_2", "Other clusters"))

fisher_surv <- survdiff(
  Surv(surv_time_years, status) ~ cluster_group,
  data = merged_data_binary
)

# Extract p-value
p_val_binary <- 1 - pchisq(fisher_surv$chisq, df = 1)
cat("\nCluster 2 vs Clusters 1/3/5/6 combined:\n")
cat(sprintf("Log-rank p-value: %.4f\n", p_val_binary))
