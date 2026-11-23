# This script calculates the inferred activity of transcription factors within each cluster based on the DEA results. 

library(dplyr)
library(tibble)
library(tidyr)
library(readr)
library(viper)
library(matrixStats)
library(decoupleR)
library(ggplot2)
library(rstatix)
library(reshape2)
library(OmnipathR)

# set working directory
set.seed(1507)
## ---- Working directory ----------------------------------------------------
## Before running, set the working directory to the "Transcriptomics/Discovery_cohort_RNAseq" folder of
## this repository (the folder that holds Inputs/ and Outputs/), e.g.
## setwd("/path/to/WD_DDLPS_phenotypes/Transcriptomics/Discovery_cohort_RNAseq")

output_dir <- "Outputs/Figures/TF_activity_violins_5C"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

metadata <- read_csv("Inputs/Metadata_discovery_WD_DD.csv")
exprs_data    <- readRDS("Inputs/Batch_corrected_data_discovery_WD_DD.rds")
clusters6 <- read.csv("Outputs/Files/Consensus_clustering/Consensus_clustering.csv", row.names = 1)
degs <- read_csv("Outputs/Files/Limma/DEGs_6cl.csv")

metadata <- metadata %>%
  left_join(clusters6, by = "sampleID")
colnames(metadata)[ncol(metadata)] <- "clusters6"

# Prepare DEG-based expression subset
degs <- as.data.frame(degs)
degs <- degs %>% distinct(gene, .keep_all = TRUE)
rownames(degs) <- degs$gene
degs <- degs[, -1]
degs <- degs[order(degs$cluster), ]

exprs_subset <- exprs_data[rownames(degs), metadata$sampleID]

# Load TF network
# CollecTRI regulon network (via OmnipathR). A saved copy in Inputs/ is used if present,
# so the analysis does not depend on the current online version of the network.
collectri_file <- "Inputs/CollecTRI_human_network.csv"
if (file.exists(collectri_file)) {
  net <- read.csv(collectri_file)
} else {
  net <- get_collectri(
    organism = "human",
    split_complexes = FALSE
  )
  write.csv(net, collectri_file, row.names = FALSE)
}

# Run VIPER
sample_acts <- run_viper(
  mat = exprs_subset,
  net = net,
  .source = "source",
  .target = "target",
  .mor = "mor",
  minsize = 5
)

# Convert to TF x samples matrix
tf_matrix <- acast(
  sample_acts,
  source ~ condition,
  value.var = "score"
)
tf_matrix <- as.data.frame(tf_matrix)
tf_matrix <- tf_matrix[, metadata$sampleID]

# Prepare TF activity table for plotting
tf_matrix_t <- t(tf_matrix) %>% as.data.frame()
tf_matrix_t$sampleID <- rownames(tf_matrix_t)
tf_matrix_t <- tf_matrix_t %>%
  left_join(metadata[, c("sampleID", "clusters6")], by = "sampleID")

# Define cluster colors
cluster_colors1 <- c(
  "C1" = "#A5C1A4",
  "C2" = "#827BB3",
  "C3" = "#DF9677",
  "C4" = "#1F70B5",
  "C5" = "#F7CBDF",
  "C6" = "#FDEDB9"
)

# Choose TF to plot
tf_of_interest <- "GLI1"

violin_df <- tf_matrix_t %>%
  select(sampleID, clusters6, all_of(tf_of_interest)) %>%
  rename(Activity = tf_of_interest) %>%
  mutate(
    Cluster = factor(
      paste0("C", clusters6),
      levels = names(cluster_colors1)
    )
  )


# Overall test: Kruskal–Wallis
kw_test <- violin_df %>%
  kruskal_test(Activity ~ Cluster)
print(kw_test)

# Wilcoxon: C4 vs all others pooled
wilcox_df <- violin_df %>%
  mutate(Group = ifelse(Cluster == "C4", "C4", "Other"))

wilcox_result <- wilcox_test(wilcox_df, Activity ~ Group) %>%
  add_significance()

print(wilcox_result)

# Add y position for bracket
y_max <- max(violin_df$Activity, na.rm = TRUE)
wilcox_result <- wilcox_result %>%
  mutate(y.position = y_max + 0.4)

# Violin plot with statistics
p <- ggplot(
  violin_df,
  aes(x = Cluster, y = Activity, fill = Cluster)
) +
  geom_violin(trim = FALSE) +
  geom_jitter(width = 0.2, size = 0.5, alpha = 0.7) +
  scale_fill_manual(values = cluster_colors1) +
  theme_minimal() +
  labs(
    title = paste0(
      "Activity of ", tf_of_interest,
      " across clusters (C4 vs rest p = ",
      signif(wilcox_result$p, 3), ")"
    ),
    x = "Cluster",
    y = "VIPER activity"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid = element_blank()
  )

print(p)

ggsave(
  filename = file.path(output_dir, paste0(tf_of_interest, "_TF_activity.pdf")),
  plot = p,
  width = 5,
  height = 4,
  units = "in")