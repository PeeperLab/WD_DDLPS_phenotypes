# Fig. 5D: copy-number OncoPrint of the discovery cohort samples with WGS (in-house) or TCGA
# copy-number data, ordered by consensus cluster.
# Fig. 5E: proportion of GLI1-amplified samples per cluster; one-sided Fisher's exact test,
# cluster 4 versus all other clusters.
set.seed(123)
## ---- Working directory ----------------------------------------------------
## Before running, set the working directory to the folder holding the
## raw WGS input data (expected folder layout: see README), e.g.
## setwd("/path/to/WGS_raw_data")

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(ggplot2)
  library(scales)
})

### input files
cna_path   <- "WGS_input/preprocessed/CNV_events.rds"
meta_file  <- "WGS_input/preprocessed/metadata_final_for_plot.tsv"
tcga_cnv   <- "WGS_input/preprocessed/data_cna.txt"
aspera_rna_meta_file <- "WGS_input/preprocessed/aspera_final_rna_metadata.csv"
clusters_k6_file     <- "clusters_k6.csv"

OUT_DIR <- "WGS_results/Figures/CNV_oncoprint_GLI1_5D_5E"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

SUPP_DIR <- file.path(OUT_DIR, "supplement_tables")
dir.create(SUPP_DIR, showWarnings = FALSE, recursive = TRUE)

## helper: robustly resolve RNA sampleID -> WGS sample_id
resolve_wgs_id <- function(x, valid_ids) {
  x <- as.character(x)

  candidates <- cbind(
    x,
    str_replace(x, "-[A-Z]\\d+$", ""),
    str_replace(x, "-R\\d+$", ""),
    str_replace(x, "-[A-Z]\\d+-R\\d+$", "")
  )

  apply(candidates, 1, function(ids) {
    hit <- ids[ids %in% valid_ids]
    if (length(hit) > 0) hit[1] else NA_character_
  })
}

read_event_table <- function(path) {
  df <- as.data.frame(readRDS(path))
  df %>%
    transmute(
      sample_id = as.character(sample_id),
      gene      = as.character(gene),
      event     = as.character(event)
    )
}

## 1) German metadata
meta <- fread(meta_file, sep = "\t", data.table = FALSE) %>%
  mutate(sample_id = as.character(sample_id)) %>%
  mutate(
    subtype = factor(as.character(subtype), ordered = TRUE),
    dataset = "German"
  )

## 1b) Load RNA metadata + consensus clusters
aspera_meta <- fread(aspera_rna_meta_file, data.table = FALSE) %>%
  transmute(
    rna_sampleID = as.character(sampleID),
    fastq_id = as.character(fastq_id)
  )

clusters_k6 <- fread(clusters_k6_file, data.table = FALSE) %>%
  transmute(
    sampleID = as.character(sampleID),
    consensus_clusters = as.character(consensus_clusters)
  )

## Only AS-* samples are used to annotate German WGS samples
cluster_as_map <- clusters_k6 %>%
  filter(str_starts(sampleID, "AS-")) %>%
  left_join(
    aspera_meta,
    by = c("sampleID" = "fastq_id")
  ) %>%
  mutate(
    wgs_sample_id = resolve_wgs_id(rna_sampleID, meta$sample_id)
  ) %>%
  filter(!is.na(wgs_sample_id)) %>%
  distinct(wgs_sample_id, .keep_all = TRUE) %>%
  transmute(
    sample_id = wgs_sample_id,
    consensus_clusters
  )

## 1c) TCGA metadata
## TCGA consensus clusters from clusters_k6.csv
tcga_cluster_map <- clusters_k6 %>%
  filter(str_starts(sampleID, "TCGA-")) %>%
  transmute(
    sample_id = sampleID,
    consensus_clusters = consensus_clusters
  ) %>%
  distinct(sample_id, .keep_all = TRUE)


## 1d) combined metadata
meta_combined <- bind_rows(cluster_as_map, tcga_cluster_map) %>%
  distinct(sample_id, .keep_all = TRUE) %>%
  filter(!is.na(consensus_clusters), consensus_clusters != "")

samples_use <- meta_combined$sample_id


## 2) Event data
## german cohort
events_german <- read_event_table(cna_path)

## TCGA cohort
tcga_data <- fread(tcga_cnv)

events_tcga <- tcga_data %>%
  rename(gene = Hugo_Symbol) %>%
  pivot_longer(
    cols = -c(gene, Entrez_Gene_Id),
    names_to = "sample_id",
    values_to = "CNV"
  ) %>%
  mutate(
    sample_id = substr(sample_id, 1, 12),
    event = case_when(
      CNV <= -2            ~ "DeepDel",
      CNV > -2 & CNV <= -1 ~ "ShallowDel",
      CNV >= 1 & CNV < 2   ~ "Gain",
      CNV >= 2             ~ "Amp",
      TRUE                 ~ NA_character_
    )
  ) %>%
  filter(!is.na(gene), gene != "", !is.na(event)) %>%
  distinct(sample_id, gene, event) %>%
  select(sample_id, gene, event)

head(events_tcga)
events_df <- bind_rows(events_german, events_tcga)


## 3) Genes list
cn_genes <- c(
  "MDM2","CDK4","HMGA2","GLI1","JUN","MAP3K5","DNM3")
genes_use <- unique(cn_genes)


## 4) long -> matrix
make_event_matrix <- function(events_df, genes, samples) {
  ev <- events_df %>%
    filter(
      gene %in% genes,
      sample_id %in% samples,
      !is.na(event), event != "",
      !is.na(gene),  gene != ""
    ) %>%
    distinct(sample_id, gene, event) %>%
    group_by(gene, sample_id) %>%
    summarise(cell_event = paste(sort(unique(event)), collapse = ";"), .groups = "drop")
  
  mat <- matrix("", nrow = length(genes), ncol = length(samples),
                dimnames = list(genes, samples))
  if (nrow(ev) == 0) return(mat)
  
  idx_g <- match(ev$gene, genes)
  idx_s <- match(ev$sample_id, samples)
  keep <- !is.na(idx_g) & !is.na(idx_s)
  mat[cbind(idx_g[keep], idx_s[keep])] <- ev$cell_event[keep]
  mat
}

pm <- make_event_matrix(events_df, genes_use, samples_use)

ensure_cols <- function(mat, cols) {
  cols <- unique(as.character(cols))
  miss <- setdiff(cols, colnames(mat))
  if (length(miss) > 0) {
    add <- matrix("", nrow = nrow(mat), ncol = length(miss),
                  dimnames = list(rownames(mat), miss))
    mat <- cbind(mat, add)
  }
  mat <- mat[, cols, drop = FALSE]
  mat
}

pm <- ensure_cols(pm, meta_combined$sample_id)


## 5) Column order + aligned annotations
cluster_map   <- setNames(as.character(meta_combined$consensus_clusters), meta_combined$sample_id)
cluster_vec <- cluster_map[colnames(pm)]

## order by consensus cluster -> sample_id
ord <- order(cluster_vec, colnames(pm), na.last = TRUE)
pm <- pm[, ord, drop = FALSE]
cluster_vec <- factor(
  cluster_map[colnames(pm)],
  levels = as.character(sort(unique(as.numeric(cluster_map[colnames(pm)]))))
)

##  Save tables 
## gene order actually used in the plot
row_order_safe <- cn_genes[cn_genes %in% rownames(pm)]
## sample annotation table (one row per plotted sample)
supp_sample_annot <- tibble(
  sample_id = colnames(pm),
  column_order = seq_along(colnames(pm)),
  consensus_clusters_raw = unname(cluster_map[colnames(pm)])
)

## long-format event table (recommended main supplementary table)
supp_event_long <- events_df %>%
  filter(
    sample_id %in% colnames(pm),
    gene %in% row_order_safe,
    !is.na(event), event != ""
  ) %>%
  distinct(sample_id, gene, event) %>%
  left_join(supp_sample_annot, by = "sample_id") %>%
  mutate(
    gene = factor(gene, levels = row_order_safe)
  ) %>%
  arrange(gene, column_order, event) %>%
  mutate(gene = as.character(gene))

## matrix table used directly for oncoprint
supp_matrix <- as.data.frame(pm, stringsAsFactors = FALSE) %>%
  rownames_to_column("gene") %>%
  mutate(gene = factor(gene, levels = row_order_safe)) %>%
  arrange(gene) %>%
  mutate(gene = as.character(gene))


## write files
fwrite(
  supp_sample_annot,
  file.path(SUPP_DIR, "Table_OncoPrint_sample_annotation.tsv"),
  sep = "\t", quote = FALSE, na = ""
)

fwrite(
  supp_event_long,
  file.path(SUPP_DIR, "Table_OncoPrint_events_long.tsv"),
  sep = "\t", quote = FALSE, na = ""
)

fwrite(
  supp_matrix,
  file.path(SUPP_DIR, "Table_OncoPrint_matrix.tsv"),
  sep = "\t", quote = FALSE, na = ""
)


## 6) Top annotation
cluster_col <- c(
  "1" = "#A5C1A4",
  "2" = "#827BB3",
  "3" = "#DF9677",
  "4" = "#1F70B5",
  "5" = "#F7CBDF",
  "6" = "#FDEDB9"
)

cluster_col <- cluster_col[names(cluster_col) %in% levels(cluster_vec)]

ha <- HeatmapAnnotation(
  Cluster = cluster_vec,
  col = list(
    Cluster = cluster_col
  ),
  annotation_name_side = "left",
  simple_anno_size = unit(7, "mm"),
  annotation_name_gp = gpar(fontsize = 14, fontface = "bold"),
   annotation_legend_param = list(
    Cluster = list(
      title_gp  = gpar(fontsize = 16, fontface = "bold"),
      labels_gp = gpar(fontsize = 16)
    )
  )
)


## 7) CNV colors & alter_fun
cnv_cols <- c(
  DeepDel    = "#044282",
  ShallowDel = "#60b1eb",
  Gain       = "#f5c66e",
  Amp        = "#E31A1C"
)

alter_fun <- list(
  background = function(x, y, w, h)
    grid.rect(x, y, w, h, gp = gpar(fill = "white", col = NA)),
  DeepDel    = function(x, y, w, h)
    grid.rect(x, y, w*0.9, h*0.85, gp = gpar(fill = cnv_cols["DeepDel"], col = NA)),
  ShallowDel = function(x, y, w, h)
    grid.rect(x, y, w*0.9, h*0.85, gp = gpar(fill = cnv_cols["ShallowDel"], col = NA)),
  Gain       = function(x, y, w, h)
    grid.rect(x, y, w*0.9, h*0.85, gp = gpar(fill = cnv_cols["Gain"], col = NA)),
  Amp        = function(x, y, w, h)
    grid.rect(x, y, w*0.9, h*0.85, gp = gpar(fill = cnv_cols["Amp"], col = NA))
)


## 8) Output
out_pdf <- file.path(OUT_DIR, "CNV_OncoPrint_with_consensus_clusters_no_subtype.pdf")
w <- max(12, ncol(pm) * 0.18)
h <- max(4, nrow(pm) * 0.2 + 2)

row_order_safe <- cn_genes[cn_genes %in% rownames(pm)]

pdf(out_pdf, width = w, height = h, useDingbats = FALSE)

ht <- oncoPrint(
  pm,
  alter_fun = alter_fun,
  col = cnv_cols,
  top_annotation = ha,
  row_order = row_order_safe,
  column_order = colnames(pm),
  column_split = data.frame(
    Cluster = cluster_vec
  ),
  column_gap = unit(1, "mm"),
  show_column_names = FALSE,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 16),
  row_names_side = "left",
  show_pct = FALSE,
  right_annotation = NULL,
  
  heatmap_legend_param = list(
    title = "CNV status",
    at = names(cnv_cols),
    labels = names(cnv_cols),
    title_gp = gpar(fontsize = 16, fontface = "bold"),
    labels_gp = gpar(fontsize = 16),
    grid_height = unit(5, "mm"),
    grid_width  = unit(5, "mm")
  ),
  
  remove_empty_rows = FALSE,
  remove_empty_columns = FALSE,
  alter_fun_is_vectorized = FALSE
)

draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

message("Saved: ", out_pdf)



####### GLI1 amplification statistics
## output file
plot_out <- file.path(OUT_DIR, "GLI1_amplified_by_cluster_stacked_barplot.pdf")

## 1) sample-level annotation table
plot_meta <- meta_combined %>%
  transmute(
    sample_id = as.character(sample_id),
    consensus_clusters = as.character(consensus_clusters)
  ) %>%
  filter(!is.na(consensus_clusters), consensus_clusters != "") %>%
  distinct(sample_id, .keep_all = TRUE)

## 2) GLI1 event table at sample level
## if a sample has Gain or Amp for GLI1 -> Amplified
## otherwise -> Not Amplified
gli1_status <- events_df %>%
  filter(gene == "GLI1") %>%
  transmute(
    sample_id = as.character(sample_id),
    event = as.character(event)
  ) %>%
  group_by(sample_id) %>%
  summarise(
    amp_status = ifelse(any(event %in% c("Gain", "Amp")), "Amplified", "Not Amplified"),
    gli1_events = paste(sort(unique(event)), collapse = ";"),
    .groups = "drop"
  )

## 3) merge back to all plotted samples
## samples without GLI1 event are treated as Not Amplified
plot_df <- plot_meta %>%
  left_join(gli1_status, by = "sample_id") %>%
  mutate(
    amp_status = ifelse(is.na(amp_status), "Not Amplified", amp_status),
    gli1_events = ifelse(is.na(gli1_events), "", gli1_events),
    consensus_clusters = factor(
      consensus_clusters,
      levels = as.character(sort(unique(as.numeric(consensus_clusters))))
    ),
    amp_status = factor(
      amp_status,
      levels = c("Not Amplified", "Amplified")
    )
  )

## 4) count + proportion per cluster
plot_stat <- plot_df %>%
  count(consensus_clusters, amp_status, name = "n") %>%
  group_by(consensus_clusters) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

## 5) Fisher test
## Cluster 4 vs all others: GLI1 amplification enrichment
test_df <- plot_df %>%
  mutate(
    cluster4_group = ifelse(as.character(consensus_clusters) == "4", "Cluster4", "Others"),
    cluster4_group = factor(cluster4_group, levels = c("Cluster4", "Others")),
    amp_status = factor(as.character(amp_status), levels = c("Amplified", "Not Amplified"))
  )

tab_2x2 <- table(test_df$cluster4_group, test_df$amp_status)

fisher_cluster4_greater <- fisher.test(tab_2x2, alternative = "greater")

p_label <- paste0("Fisher p = ", signif(fisher_cluster4_greater$p.value, 3))

## 6) stacked barplot
p <- ggplot(plot_stat, aes(x = consensus_clusters, y = prop, fill = amp_status)) +
  geom_col(width = 0.9, color = NA) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1.10),
    expand = c(0, 0)
  ) +
  scale_fill_manual(
    values = c(
      "Not Amplified" = "#BDBDBD",
      "Amplified" = "#FF2D0A"
    ),
    name = "Amplification"
  ) +
  labs(
    title = "Proportion of GLI1-amplified Samples per Cluster",
    x = "Cluster",
    y = "Proportion of samples"
  ) +
  annotate(
    "text",
    x = mean(seq_along(levels(plot_df$consensus_clusters))),
    y = 1.04,
    label = p_label,
    size = 7
  ) +
  theme_classic(base_size = 18) +
  theme(
    plot.title = element_text(size = 26, face = "plain"),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 18),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18)
  )

print(p)

ggsave(plot_out, p, width = 12, height = 8, units = "in")

supp_gli1 <- plot_df %>%
  arrange(consensus_clusters, sample_id)

fwrite(
  supp_gli1,
  file.path(OUT_DIR, "Table_GLI1_amplification_by_cluster.tsv"),
  sep = "\t",
  quote = FALSE,
  na = ""
)


prop_summary <- test_df %>%
  count(cluster4_group, amp_status) %>%
  tidyr::pivot_wider(names_from = amp_status, values_from = n, values_fill = 0) %>%
  mutate(
    total = Amplified + `Not Amplified`,
    amp_prop = Amplified / total
  )

stat_out <- tibble::tibble(
  comparison = "Cluster4_vs_Others",
  amplified_cluster4 = tab_2x2["Cluster4", "Amplified"],
  not_amplified_cluster4 = tab_2x2["Cluster4", "Not Amplified"],
  amplified_others = tab_2x2["Others", "Amplified"],
  not_amplified_others = tab_2x2["Others", "Not Amplified"],
  p_value_one_sided = fisher_cluster4_greater$p.value,
  odds_ratio = unname(fisher_cluster4_greater$estimate),
  conf_low = fisher_cluster4_greater$conf.int[1],
  conf_high = fisher_cluster4_greater$conf.int[2]
)

fwrite(
  as.data.frame.matrix(tab_2x2, stringsAsFactors = FALSE) %>%
    tibble::rownames_to_column("group"),
  file.path(OUT_DIR, "Table_GLI1_cluster4_vs_others_2x2.tsv"),
  sep = "\t",
  quote = FALSE,
  na = ""
)

fwrite(
  prop_summary,
  file.path(OUT_DIR, "Table_GLI1_cluster4_vs_others_proportions.tsv"),
  sep = "\t",
  quote = FALSE,
  na = ""
)

fwrite(
  stat_out,
  file.path(OUT_DIR, "Table_GLI1_cluster4_vs_others_fisher.tsv"),
  sep = "\t",
  quote = FALSE,
  na = ""
)