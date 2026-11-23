
#!/usr/bin/env Rscript
# Fig. 1F: OncoPrint of the in-house WGS samples (hg19): selected mutations, copy-number
# events and structural variants, with TMB and FGA, split by histological subtype.
# Inputs: outputs of S4 (MUT_SV_events.rds) and S5 (CNV_events.rds, metadata_final_for_plot.tsv).
set.seed(123)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

## Working directory and inputs
## ---- Working directory ----------------------------------------------------
## Before running, set the working directory to the folder holding the
## raw WGS input data (expected folder layout: see README), e.g.
## setwd("/path/to/WGS_raw_data")

mut_sv_path <- "WGS_input/preprocessed/MUT_SV_events.rds"
cna_path    <- "WGS_input/preprocessed/CNV_events.rds"
meta_file   <- "WGS_input/preprocessed/metadata_final_for_plot.tsv"

OUT_DIR <- "WGS_results/Figures/OncoPrint_1F"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

## Read inputs
events_mut_sv <- readRDS(mut_sv_path)
events_cna    <- readRDS(cna_path)
meta <- fread(meta_file, sep = "\t", data.table = FALSE) %>% as.data.frame()
meta$sample_id <- as.character(meta$sample_id)

desired_levels <- c("DDLPS","WD/DDLPS","WDLPS","PLPS","MLPS")
meta <- meta %>%
  mutate(
    subtype = factor(as.character(subtype), levels = desired_levels, ordered = TRUE),
    patient_id = as.character(patient_id)
  ) %>%
  arrange(subtype, patient_id, sample_id)

samples_use <- meta$sample_id

# Genes to show 
mut_genes <- c("TP53","ATRX","RB1","PTEN")
cn_genes  <- c("MDM2","CDK4","HMGA2","GLI1","TERT","DNM3","JUN","MAP3K5","YAP1","CDKN2A","PTEN","TP53","RB1")
sv_genes  <- c("DDIT3", "FUS", "CSMD1")

genes_use <- unique(c(mut_genes, cn_genes, sv_genes))

###Build OncoPrint matrix (gene x sample), each cell has "event1;event2" 
events_all <- bind_rows(events_mut_sv, events_cna) %>%
  transmute(
    sample_id = as.character(sample_id),
    gene      = as.character(gene),
    event     = as.character(event)
  ) %>%
  filter(gene %in% genes_use, !is.na(event), event != "", !is.na(gene), gene != "") %>%
  distinct(sample_id, gene, event)

cell_events <- events_all %>%
  group_by(sample_id, gene) %>%
  summarise(cell_event = paste(sort(unique(event)), collapse = ";"), .groups = "drop")

mat <- matrix("", nrow = length(genes_use), ncol = length(samples_use),
              dimnames = list(genes_use, samples_use))

for (k in seq_len(nrow(cell_events))) {
  g <- cell_events$gene[k]
  s <- cell_events$sample_id[k]
  e <- cell_events$cell_event[k]
  if (g %in% rownames(mat) && s %in% colnames(mat)) mat[g, s] <- e
}

# Top annotations vectors (aligned to mat columns) 
subtype_vec <- setNames(meta$subtype, meta$sample_id)[colnames(mat)]

get_num_vec <- function(df, col, samples) {
  if (!col %in% colnames(df)) return(setNames(rep(NA_real_, length(samples)), samples))
  v <- suppressWarnings(as.numeric(df[[col]]))
  setNames(v, df$sample_id)[samples]
}

tcc_vec    <- get_num_vec(meta, "tcc",    colnames(mat))
ploidy_vec <- get_num_vec(meta, "ploidy", colnames(mat))
fga_vec    <- get_num_vec(meta, "FGA",    colnames(mat))
tmb_vec    <- get_num_vec(meta, "TMB",    colnames(mat))

## Colors 
safe_range <- function(x, pad = 0.05) {
  r <- range(x, na.rm = TRUE)
  if (!any(is.finite(r))) return(c(0, 1))
  if (r[1] == r[2]) r <- r + c(-pad, pad)
  r
}

legend_breaks <- function(x, n = 4) {
  r <- range(x, na.rm = TRUE)
  if (!any(is.finite(r))) return(NULL)
  pretty(r, n = n)
}


sub_col <- c(DDLPS = "#F18072",
 WDLPS = "#80B1D2",
 PLPS = "#060771",
 MLPS = "#4F7942",
 `WD/DDLPS` = "#F7CBDF")

# Order columns by subtype
ord <- order(subtype_vec, na.last = TRUE)
mat <- mat[, ord, drop = FALSE]

subtype_vec <- subtype_vec[ord]
tmb_vec     <- tmb_vec[ord]
fga_vec     <- fga_vec[ord]

col_order <- colnames(mat)

# Top annotation 
ha <- HeatmapAnnotation(
  Subtype = subtype_vec,
  TMB     = anno_barplot(tmb_vec, border = FALSE, 
                         axis_param = list(gp = gpar(fontsize = 11))),
  FGA     = anno_barplot(fga_vec, border = FALSE, 
                         axis_param = list(gp = gpar(fontsize = 11))),
  col = list(Subtype = sub_col),
  annotation_name_side = "left",
  
  annotation_name_gp = gpar(fontsize = 14,
                            fontface = "bold"
  ),
  
  annotation_legend_param = list(
    Subtype = list(
      title_gp = gpar(fontsize = 16,
                      fontface = "bold"
      ), 
      labels_gp = gpar(fontsize = 16)
    )
  )
)

## alter_fun & color mapping 
alter_fun <- list(
  background = function(x, y, w, h) grid.rect(x, y, w, h, gp=gpar(fill="white", col=NA)),
  
  Missense   = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.3, gp=gpar(fill="#2CA25F", col=NA)),
  Nonsense   = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.3, gp=gpar(fill="#6ef5ec", col=NA)),
  Frameshift = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.3, gp=gpar(fill="#af6ef5", col=NA)),
  Inframe    = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.3, gp=gpar(fill="#A50F15", col=NA)),
  Splice     = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.3, gp=gpar(fill="#6A3D9A", col=NA)),
  SNV        = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.3, gp=gpar(fill="#2CA25F", col=NA)),
  INDEL      = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.3, gp=gpar(fill="#FB9A29", col=NA)),
  
  DeepDel    = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp=gpar(fill="#044282", col=NA)),
  ShallowDel = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp=gpar(fill="#60b1eb", col=NA)),
  Gain       = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp=gpar(fill="#f5c66e", col=NA)),
  Amp        = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp=gpar(fill="#E31A1C", col=NA)),
  
  SV_DEL     = function(x, y, w, h) { dx <- w*0.45; dy1 <- h*0.30; dy2 <- h*0.35
  grid.polygon(unit.c(x - dx, x + dx, x),
               unit.c(y - dy1, y - dy1, y + dy2), gp=gpar(fill="#FF7F00", col=NA))
  },
  SV_DUP     = function(x, y, w, h) { dx <- w*0.45; dy1 <- h*0.30; dy2 <- h*0.35
  grid.polygon(unit.c(x - dx, x + dx, x),
               unit.c(y - dy1, y - dy1, y + dy2), gp=gpar(fill="#1F78B4", col=NA))
  },
  SV_INV     = function(x, y, w, h) { dx <- w*0.45; dy1 <- h*0.30; dy2 <- h*0.35
  grid.polygon(unit.c(x - dx, x + dx, x),
               unit.c(y - dy1, y - dy1, y + dy2), gp=gpar(fill="#ec6ef5", col=NA))
  },
  SV_TRA     = function(x, y, w, h) { dx <- w*0.45; dy1 <- h*0.30; dy2 <- h*0.35
  grid.polygon(unit.c(x - dx, x + dx, x),
               unit.c(y - dy1, y - dy1, y + dy2), gp=gpar(fill="#33A02C", col=NA))
  },
  SV         = function(x, y, w, h) { dx <- w*0.45; dy1 <- h*0.30; dy2 <- h*0.35
  grid.polygon(unit.c(x - dx, x + dx, x),
               unit.c(y - dy1, y - dy1, y + dy2), gp=gpar(fill="#E31A1C", col=NA))
  }
)

mutation_cols <- c(Missense="#2CA25F", Nonsense="#6ef5ec", Frameshift="#af6ef5", Inframe="#A50F15", Splice="#6A3D9A", SNV="#2CA25F", INDEL="#FB9A29")
cnv_cols      <- c(DeepDel="#044282", ShallowDel="#60b1eb", Gain="#f5c66e", Amp="#E31A1C", LOH="#636363")
sv_cols       <- c(SV_DEL="#FF7F00", SV_DUP="#1F78B4", SV_INV="#ec6ef5", SV_TRA="#33A02C", SV="#E31A1C")

split_group <- function(m, ev_set, genes_keep) {
  genes_keep <- intersect(unique(genes_keep), rownames(m))
  if (length(genes_keep) == 0) return(NULL)
  
  m_sub <- m[genes_keep, , drop = FALSE]
  out <- matrix("", nrow=nrow(m_sub), ncol=ncol(m_sub), dimnames=dimnames(m_sub))
  
  for (i in seq_len(nrow(m_sub))) {
    for (j in seq_len(ncol(m_sub))) {
      hits <- unlist(strsplit(m_sub[i, j], ";", fixed = TRUE))
      hits <- hits[hits %in% ev_set]
      out[i, j] <- paste(hits, collapse = ";")
    }
  }
  if (all(out == "")) return(NULL)
  out
}

panel_defs <- list(
  list(label="Mutation", events=names(mutation_cols), cols=mutation_cols, legend="Mutation", genes=mut_genes),
  list(label="CNV", events=names(cnv_cols),      cols=cnv_cols,      legend="CNV",      genes=cn_genes),
  list(label="SV",  events=names(sv_cols),       cols=sv_cols,       legend="SV",       genes=sv_genes)
)

## Output
out_pdf <- file.path(OUT_DIR, "DNA_landscape_OncoPrint_hg19.pdf")
w <- max(12, ncol(mat) * 0.18)
h <- max(8, nrow(mat) * 0.25 + 2)

pdf(out_pdf, width = w, height = h, useDingbats = FALSE)

ht_list <- NULL
first_panel <- TRUE
gap_id <- 1

for (pd in panel_defs) {
  pm <- split_group(mat, pd$events, pd$genes)
  if (is.null(pm)) next
  
  # panel-specific row order (stable)
  ro <- match(pd$genes, rownames(pm))
  ro <- ro[!is.na(ro)]
  
  ht <- oncoPrint(
    pm,
    alter_fun = alter_fun,
    col = pd$cols,
    top_annotation = if (first_panel) ha else NULL,
    row_order = ro,
    column_order = colnames(pm),      
    column_split = subtype_vec,       # split by subtype
    column_gap   = unit(1, "mm"),
    # border = "#E3E3E3",              
    # border_gp = gpar(lwd = 1),  
    column_title_gp = gpar(
      #fontface="bold", 
      fontsize = 14),
    row_title = pd$label,
    row_title_side = "left",
    row_title_rot = 90,
    row_title_gp = gpar(
      fontface="bold", 
      fontsize=16),
    show_column_names = FALSE,
    show_row_names = TRUE,
    show_pct = FALSE,
    right_annotation = NULL,  ## this will remove the barplot on the right
    row_names_gp = gpar(fontsize=16),
    heatmap_legend_param = list(
      title = pd$legend,
      at = names(pd$cols),       
      labels = names(pd$cols),
      title_gp = gpar(fontsize = 16, fontface = "bold"),
      labels_gp = gpar(fontsize = 16),
      grid_height = unit(6, "mm"), 
      grid_width = unit(6, "mm")
    ),
    remove_empty_rows = TRUE,
    remove_empty_columns = FALSE,
    alter_fun_is_vectorized = FALSE
  )
  
  if (is.null(ht_list)) {
    ht_list <- ht
  } else {
    gap_ht <- Heatmap(matrix(0, nrow=1, ncol=ncol(pm)),
                      name=paste0("gap", gap_id), col=c("0"= "#E3E3E3"),
                      show_heatmap_legend=FALSE, show_row_names=FALSE, show_column_names=FALSE,
                      height=unit(0.3, "mm"), border=NA)
    ht_list <- ht_list %v% gap_ht %v% ht
    gap_id <- gap_id + 1
  }
  first_panel <- FALSE
}

if (is.null(ht_list)) stop("No oncoPrint panels to draw (empty matrices).")

draw(ht_list,
     heatmap_legend_side = "right",
     annotation_legend_side = "right"
     # annotation_legend_list = ann_legends
)

dev.off()
message("Saved: ", out_pdf)



