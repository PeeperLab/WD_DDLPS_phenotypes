set.seed(123)
## ---- Working directory ----------------------------------------------------
## Before running, set the working directory to the folder holding the
## raw WGS input data (expected folder layout: see README), e.g.
## setwd("/path/to/WGS_raw_data")
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(GenomicRanges)
  library(IRanges)
  library(rtracklayer)
  library(tidyr)
})

### input
cnv_dir        <- "WGS_input/renamed/common/cnv"
gtf_path       <- "WGS_input/gencode.v19.annotation.gtf"
meta_file      <- "WGS_input/preprocessed/metadata_with_TMB.tsv"
manifest_file  <- "WGS_input/preprocessed/WGS_data_manifest.tsv"

OUT_DIR <- "WGS_input/preprocessed"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
### helper
std_chr <- c(as.character(1:22), "X", "Y")

norm_chr <- function(x) {
  x <- gsub("^chr", "", x, ignore.case = TRUE)
  toupper(as.character(x))
}


keep_std_chr_gr <- function(gr) {
  gr[as.character(seqnames(gr)) %in% std_chr]
}

load_genes_from_gtf <- function(gtf_path) {
  g <- rtracklayer::import(gtf_path)
  g <- g[g$type == "gene"]

  gene_name <- mcols(g)$gene_name
  gene_id   <- mcols(g)$gene_id
  if (is.null(gene_name)) gene_name <- gene_id
  if (is.null(gene_id)) gene_id <- gene_name

  chr_vec <- norm_chr(as.character(seqnames(g)))
  keep <- !is.na(chr_vec) & chr_vec %in% std_chr

  out <- GRanges(
    seqnames  = chr_vec[keep],
    ranges    = IRanges(start = start(g)[keep], end = end(g)[keep]),
    gene_name = as.character(gene_name)[keep],
    gene_id   = as.character(gene_id)[keep]
  )
  seqlevels(out) <- std_chr
  out
}

manifest <- fread(manifest_file, sep = "\t", data.table = FALSE)

####### read and process CNV data
#### ACE CNV reader + CNA burden + WGD
read_ace_most_important <- function(path, drop_nocall = TRUE, eps = 1e-8) {
  hdr <- readLines(path, n = 50, warn = FALSE)

  get_kv_num <- function(key) {
    line <- hdr[grepl(paste0("^#", key, ":"), hdr)]
    as.numeric(sub(paste0("^#", key, ":"), "", line[1]))
  }

  tcc         <- get_kv_num("tcc")
  ploidy      <- get_kv_num("ploidy")
  roundPloidy <- get_kv_num("roundPloidy")
  quality     <- get_kv_num("quality")

  seg <- fread(
    path,
    sep = "\t",
    header = TRUE,
    skip = "#chromosome",
    data.table = FALSE,
    showProgress = FALSE
  )

  names(seg)[names(seg) == "#chromosome"] <- "chromosome"

  seg$chromosome    <- norm_chr(seg$chromosome)
  seg$start         <- as.integer(seg$start)
  seg$end           <- as.integer(seg$end)
  seg$length        <- as.numeric(seg$length)
  seg$covRatio      <- as.numeric(seg$covRatio)
  seg$TCN           <- as.numeric(seg$TCN)
  seg$c1Mean        <- as.numeric(seg$c1Mean)
  seg$c2Mean        <- as.numeric(seg$c2Mean)
  seg$CNA.type      <- as.character(seg$CNA.type)
  seg$NbrOfHetsSNPs <- as.integer(seg$NbrOfHetsSNPs)

  seg <- seg[
    seg$chromosome %in% std_chr &
      !is.na(seg$start) & !is.na(seg$end) &
      seg$start > 0 & seg$end > 0 & seg$start <= seg$end &
      is.finite(seg$TCN),
    , drop = FALSE
  ]

  seg$is_nocall <-
    (is.finite(seg$covRatio) & abs(seg$covRatio) <= eps) &
    (is.finite(seg$TCN) & abs(seg$TCN) <= eps) &
    (is.na(seg$NbrOfHetsSNPs) | seg$NbrOfHetsSNPs <= 0L)

  n_nocall <- sum(seg$is_nocall, na.rm = TRUE)
  bp_nocall <- sum(seg$length[seg$is_nocall], na.rm = TRUE)

  if (isTRUE(drop_nocall)) {
    seg <- seg[!seg$is_nocall, , drop = FALSE]
  }

  list(
    header = data.frame(
      tcc = tcc,
      ploidy = ploidy,
      roundPloidy = roundPloidy,
      quality = quality,
      stringsAsFactors = FALSE
    ),
    seg = seg,
    qc = data.frame(n_nocall = n_nocall, bp_nocall = bp_nocall, stringsAsFactors = FALSE)
  )
}

### calculate read depth/coverage of baseline
infer_covratio_baseline <- function(seg, min_len = 1e6) {
  v <- seg$covRatio
  okv <- is.finite(v) & v > 0 & is.finite(seg$length) & seg$length > 0
  if (!any(okv)) return(NA_real_)

  ct <- tolower(trimws(seg$CNA.type))
  is_neu <- ct %in% c("tcnneutral", "neutral", "tcn_neutral", "tcn-neutral")

  cand <- v[okv & is_neu & seg$length >= min_len]
  if (length(cand) >= 10) return(median(cand))

  cand2 <- v[okv & is_neu]
  if (length(cand2) >= 10) return(median(cand2))

  median(v[okv])
}
### calculate seg mean from covratio, log2(seg/baseline), cap level: -10 and 10
compute_segmean_from_covratio <- function(seg, clip_low = -10, clip_high = 10) {
  base <- infer_covratio_baseline(seg)
  if (!is.finite(base) || base <= 0) return(rep(NA_real_, nrow(seg)))

  x <- seg$covRatio
  y <- rep(NA_real_, nrow(seg))

  ok <- is.finite(x) & x > 0
  y[ok] <- log2(x[ok] / base)
  y[ok] <- pmin(pmax(y[ok], clip_low), clip_high)
  y
}
## calculate fraction of genome altered (FGA)
## converts the log2 ratio of CNV segments into the proportion of the genome covered by gain/deletion/any alteration (FGA, GainFrac, LossFrac)
compute_cna_burden_from_segmean <- function(seg, segmean, thr = 0.2) {
  ok <- is.finite(seg$length) & seg$length > 0 & is.finite(segmean)
  if (!any(ok)) return(list(FGA = NA_real_, GainFrac = NA_real_, LossFrac = NA_real_))

  len <- seg$length[ok]
  total <- sum(len)
  if (total <= 0) return(list(FGA = NA_real_, GainFrac = NA_real_, LossFrac = NA_real_))

  list(
    FGA = sum(len[abs(segmean[ok]) > thr]) / total,
    GainFrac = sum(len[segmean[ok] > thr]) / total,
    LossFrac = sum(len[segmean[ok] < -thr]) / total
  )
}

#### 4-level CNV calling
call_cna_event_4level <- function(CNA.type, TCN, Ploidy,
                                  deep_del_rel = 0.25,
                                  shallow_del_rel = 0.75,
                                  gain_rel = 1.25,
                                  amp_rel = 2) {
  TCN <- as.numeric(TCN)
  Ploidy <- as.numeric(Ploidy)
  if (!is.finite(TCN) || !is.finite(Ploidy) || Ploidy <= 0) return(NA_character_)

  rel <- TCN / Ploidy
  ct  <- tolower(as.character(CNA.type))
  if (str_detect(ct, "homodel") || TCN <= 0.1 || rel <= deep_del_rel) return("DeepDel")
  if (rel < shallow_del_rel) return("ShallowDel")
  if (rel >= amp_rel) return("Amp")
  if (rel >= gain_rel) return("Gain")
  NA_character_
}

weighted_mean <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  sum(x[ok] * w[ok]) / sum(w[ok])
}
# Gene-level CNA from seg/hdr 
cna_events_from_seg_weighted <- function(seg, hdr, genes_gr_sub,
                                         loh_min_frac = 0.5,
                                         hdel_min_frac = 0.5) {
  P <- hdr$ploidy[1]
  if (!is.finite(P) || P <= 0) {
    P <- hdr$roundPloidy[1]
    if (!is.finite(P) || P <= 0) P <- 2
  }

  seg_gr <- GRanges(
    seqnames = seg$chromosome,
    ranges = IRanges(start = as.integer(seg$start), end = as.integer(seg$end)),
    TCN = as.numeric(seg$TCN),
    CNA.type = as.character(seg$CNA.type)
  )
  seg_gr <- keep_std_chr_gr(seg_gr)

  hits <- findOverlaps(genes_gr_sub, seg_gr, ignore.strand = TRUE)
  if (length(hits) == 0) {
    return(list(events = data.frame(), gene_level = data.frame(), ploidy_used = P))
  }

  ov <- pintersect(ranges(genes_gr_sub[queryHits(hits)]), ranges(seg_gr[subjectHits(hits)]))
  ovw <- width(ov)

  df0 <- data.frame(
    gene = mcols(genes_gr_sub)$gene_name[queryHits(hits)],
    TCN = mcols(seg_gr)$TCN[subjectHits(hits)],
    CNA.type = mcols(seg_gr)$CNA.type[subjectHits(hits)],
    ovw = ovw,
    stringsAsFactors = FALSE
  ) %>%
    filter(is.finite(ovw) & ovw > 0, !is.na(gene), gene != "")

  if (nrow(df0) == 0) {
    return(list(events = data.frame(), gene_level = data.frame(), ploidy_used = P))
  }

  gene_level <- df0 %>%
    mutate(
      ct = tolower(trimws(CNA.type)),
      is_hdel = str_detect(ct, "homodel")
    ) %>%
    group_by(gene) %>%
    summarise(
      ovw_sum = sum(ovw, na.rm = TRUE),
      TCN_w = weighted_mean(TCN, ovw),
      hdel_frac = ifelse(ovw_sum > 0, sum(ovw[is_hdel], na.rm = TRUE) / ovw_sum, NA_real_),
      .groups = "drop"
    )

  gene_events <- gene_level %>%
    mutate(
      cn_event = mapply(
        call_cna_event_4level,
        CNA.type = ifelse(is.finite(hdel_frac) & hdel_frac >= hdel_min_frac, "HomoDel", ""),
        TCN = TCN_w,
        MoreArgs = list(Ploidy = P),
        USE.NAMES = FALSE
      )
    ) %>%
    select(gene, cn_event) %>%
    tidyr::pivot_longer(cols = c(cn_event), names_to = "kind", values_to = "event") %>%
    filter(!is.na(event), event != "") %>%
    distinct(gene, event)

  list(events = gene_events, gene_level = gene_level, ploidy_used = P)
}

genes_gr <- load_genes_from_gtf(gtf_path)

cnv_res <- lapply(seq_len(nrow(manifest)), function(i) {
  sid <- as.character(manifest$sample_id[i])

  obj <- read_ace_most_important(manifest$cnv_path[i])
  seg <- obj$seg
  hdr <- obj$header

  segmean <- compute_segmean_from_covratio(seg)
  burden <- compute_cna_burden_from_segmean(seg, segmean, thr = 0.2)
  # Use all genes from the GTF for CNV events
  ev_obj <- cna_events_from_seg_weighted(seg, hdr, genes_gr)
  ev <- ev_obj$events %>% mutate(sample_id = sid, datatype = "CNA")

  hd <- hdr %>%
    mutate(
      sample_id = sid,
      FGA = burden$FGA,
      GainFrac = burden$GainFrac,
      LossFrac = burden$LossFrac
    )

  list(events = ev, header = hd)
})

events_cna <- bind_rows(lapply(cnv_res, `[[`, "events")) %>%
  transmute(
    sample_id = as.character(sample_id),
    gene      = as.character(gene),
    event     = as.character(event),
    datatype  = as.character(datatype)
  )

cnv_hdr_df <- bind_rows(lapply(cnv_res, `[[`, "header"))

saveRDS(events_cna, file.path(OUT_DIR, "CNV_events.rds"))
fwrite(events_cna, file.path(OUT_DIR, "CNV_events.tsv"), sep = "\t")

# saveRDS(cnv_hdr_df, file.path(OUT_DIR, "CNV_headr.rds"))
# fwrite(cnv_hdr_df, file.path(OUT_DIR, "CNV_headr.tsv"), sep = "\t")

metadata <- fread(meta_file, sep = "\t", data.table = FALSE)
metadata$sample_id <- as.character(metadata$sample_id)
metadata <- metadata[, !colnames(metadata) %in% c("tcc", "ploidy")]

metadata_final <- metadata %>%
  left_join(cnv_hdr_df, by = "sample_id")

fwrite(metadata_final, file.path(OUT_DIR, "metadata_final_for_plot.tsv"), sep = "\t")
