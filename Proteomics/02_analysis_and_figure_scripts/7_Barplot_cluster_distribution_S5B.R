# This script visualizes the distribution of clusters of multiple samples from the same tumor. 
# Only tumors with at least 3 samples are included. 

library(openxlsx)  # read.xlsx()
library(dplyr)
library(ggplot2)

set.seed(1507)
## ---- Working directory ----------------------------------------------------
## Before running, set the working directory to the "Proteomics" folder of
## this repository (the folder that holds Inputs/ and Outputs/), e.g.
## setwd("/path/to/WD_DDLPS_phenotypes/Proteomics")

# Load metadata and cluster info
proteomic_metadata <- read.xlsx("Inputs/Proteomics_metadata.xlsx")
clusters <- read.csv("Outputs/Files/Consensus_clustering/Consensus_clustering.csv", row.names = 1)     

output_dir <- "Outputs/Figures/Cluster_barplots_S5B"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Put ID as rownames
rownames(proteomic_metadata) <- proteomic_metadata[[1]]
clusters$ID <- rownames(clusters)

# Merge cluster info into metadata
df <- proteomic_metadata %>%
  select(ID, pat) %>%
  left_join(clusters[, c("ID", "consensus_clusters")], by = "ID") %>%
  rename(final_cluster = consensus_clusters)

df <- df %>%
  filter(!is.na(final_cluster) & !final_cluster %in% c("Healthy_Fat", "Lipoma"))

# Filter patients with at least 3 samples
df_filtered <- df %>%
  group_by(pat) %>%
  filter(n() >= 3) %>%
  ungroup()

# Compute counts and proportions per patient and cluster
df_prop <- df_filtered %>%
  group_by(pat, final_cluster) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(pat) %>%
  mutate(prop = count / sum(count)) %>%
  ungroup()

# Treat clusters as factor and define custom order/colors
cluster_levels <- c("1","2","3","4","5","6")
cluster_colors <- c(
  "1" = "#A5C1A4",
  "2" = "#827BB3",
  "3" = "#E34363",
  "4" = "#DF9677",
  "5" = "#1F70B5",
  "6" = "#FDEDB9"
)

df_prop$final_cluster <- factor(df_prop$final_cluster, levels = cluster_levels)

# Define cluster order
cluster_priority <- c("6","1","3","5","4","2")

# Compute number of clusters and top cluster per patient
patient_priority <- df_filtered %>%
  group_by(pat) %>%
  summarise(
    clusters_present = list(unique(final_cluster)),
    n_clusters = n_distinct(final_cluster),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    # Find the earliest cluster in priority that is present
    top_cluster = cluster_priority[min(which(cluster_priority %in% clusters_present))]
  ) %>%
  ungroup() %>%
  arrange(n_clusters, match(top_cluster, cluster_priority))  # first by n_clusters, then by cluster priority

# Extract patient order
patient_order <- patient_priority$pat

# Apply to df_prop
df_prop$pat <- factor(df_prop$pat, levels = patient_order)

# Plot
p <- ggplot(df_prop, aes(x = pat, y = prop, fill = final_cluster)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = cluster_colors) +
  labs(
    x = "Patient",
    y = "Proportion of Samples",
    fill = "Cluster",
    title = "Proportion of Clusters per Patient with ≥ 3 samples"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

print(p)

ggsave(
  filename = file.path(output_dir,"Cluster_barplot.pdf"),
  plot = p,
  width = 6,
  height = 4,
  units = "in",
)

