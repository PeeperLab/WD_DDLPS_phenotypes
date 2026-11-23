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
})

HAS_MUTPAT <- requireNamespace("MutationalPatterns", quietly = TRUE)
HAS_BSGEN  <- requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)
if (HAS_MUTPAT) library(MutationalPatterns)
if (HAS_BSGEN) library(BSgenome.Hsapiens.UCSC.hg19)

manifest_file <- "WGS_input/preprocessed/WGS_data_manifest.tsv"
meta_file     <- "WGS_input/preprocessed/dna_metadata_ordered.tsv"
gtf_path      <- "WGS_input/gencode.v19.annotation.gtf"
OUT_DIR       <- "WGS_input/preprocessed"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

## helper
`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b

std_chr <- c(as.character(1:22), "X", "Y")
norm_chr <- function(x) {
  x <- gsub("^chr", "", x, ignore.case = TRUE)
  toupper(as.character(x))
}

to_ucsc_chr <- function(x) {
  paste0("chr", norm_chr(x))
}

### manifest for file path
manifest <- fread(manifest_file, sep = "\t", data.table = FALSE)

classify_mut <- function(func, exo, ref, alt) {
  func <- tolower(func %||% "")
  exo  <- tolower(exo %||% "")

  if (str_detect(exo, "synonymous")) return(NA_character_)
  if (str_detect(func, "splic") | str_detect(exo, "splice")) return("Splice")
  if (str_detect(exo, "frameshift")) return("Frameshift")
  if (str_detect(exo, "stopgain|stop_gained|nonsense")) return("Nonsense")
  if (str_detect(exo, "nonframeshift|inframe")) return("Inframe")
  if (str_detect(exo, "nonsynonymous|missense")) return("Missense")

  if (!is.na(ref) && !is.na(alt) && nchar(ref) == 1 && nchar(alt) == 1) return("SNV")
  "INDEL"
}

### read the events in SNV and Indel table
read_snv_indel_table_events <- function(path, sample_id) {
  dt <- tryCatch(
    fread(path, sep = "\t", data.table = FALSE, quote = ""),
    error = function(e) {
      warning("SNV/INDEL read failed: ", basename(path), " | ", conditionMessage(e))
      return(NULL)
    }
  )
  if (is.null(dt) || nrow(dt) == 0) return(data.frame())

  req <- c("FILTER", "ANNOVAR_FUNCTION", "EXONIC_CLASSIFICATION", "REF", "ALT", "GENE")
  miss <- setdiff(req, names(dt))
  if (length(miss) > 0) {
    warning("SNV/INDEL skip missing columns: ", basename(path), " | ", paste(miss, collapse = ","))
    return(data.frame())
  }

  dt <- dt[is.na(dt$FILTER) | dt$FILTER %in% c("PASS", ".", ""), , drop = FALSE]
  if (nrow(dt) == 0) return(data.frame())

  func <- as.character(dt$ANNOVAR_FUNCTION)
  exo  <- as.character(dt$EXONIC_CLASSIFICATION)
  ref  <- as.character(dt$REF)
  alt  <- as.character(dt$ALT)

  keep <- str_detect(tolower(func), "exonic|splic") |
    str_detect(tolower(exo), "splice|frameshift|stop|missense|nonsynonymous|inframe|nonframeshift|synonymous")

  dt   <- dt[keep, , drop = FALSE]
  func <- func[keep]
  exo  <- exo[keep]
  ref  <- ref[keep]
  alt  <- alt[keep]
  if (nrow(dt) == 0) return(data.frame())

  event <- mapply(classify_mut, func, exo, ref, alt, USE.NAMES = FALSE)
  gene_raw <- as.character(dt$GENE)

  ok <- !is.na(event) & event != "" &
    !is.na(gene_raw) & gene_raw != "" & gene_raw != "." & gene_raw != "NA"
  gene_raw <- gene_raw[ok]
  event    <- event[ok]
  if (length(gene_raw) == 0) return(data.frame())

  gene_list <- strsplit(gene_raw, "[,;|&]")
  n_each <- vapply(gene_list, length, integer(1))
  gene_vec  <- trimws(unlist(gene_list, use.names = FALSE))
  event_vec <- rep(event, times = n_each)

  out <- data.frame(
    sample_id = rep(sample_id, length(gene_vec)),
    gene      = gene_vec,
    event     = event_vec,
    datatype  = "MUT",
    stringsAsFactors = FALSE
  )

  out %>%
    filter(!is.na(gene), gene != "", gene != ".", !is.na(event), event != "") %>%
    distinct(sample_id, gene, event, datatype)
}

read_sv_events <- function(path, sample_id, min_eventScore = 4, min_minEvidence = 5) {
  dt <- tryCatch(
    fread(path, sep = "\t", data.table = FALSE, quote = ""),
    error = function(e) {
      warning("SV read failed: ", basename(path), " | ", conditionMessage(e))
      return(NULL)
    }
  )
  if (is.null(dt) || nrow(dt) == 0) return(data.frame())

  req <- c("svtype", "eventScore", "evidence1", "evidence2", "gene1", "gene2")
  miss <- setdiff(req, names(dt))
  if (length(miss) > 0) {
    warning("SV skip missing columns: ", basename(path), " | ", paste(miss, collapse = ","))
    return(data.frame())
  }

  svt <- toupper(as.character(dt$svtype))
  event <- dplyr::case_when(
    stringr::str_detect(svt, "DEL") ~ "SV_DEL",
    stringr::str_detect(svt, "DUP") ~ "SV_DUP",
    stringr::str_detect(svt, "INV") ~ "SV_INV",
    stringr::str_detect(svt, "TRA|BND") ~ "SV_TRA",
    TRUE ~ "SV"
  )

  evScore <- suppressWarnings(as.numeric(dt$eventScore))
  e1 <- suppressWarnings(as.numeric(dt$evidence1))
  e2 <- suppressWarnings(as.numeric(dt$evidence2))
  minEv <- pmin(e1, e2)

  keep <- which(is.finite(evScore) & evScore >= min_eventScore &
                  is.finite(minEv) & minEv >= min_minEvidence)
  if (length(keep) == 0) return(data.frame())

  parse_gene_field <- function(x) {
    entries <- trimws(unlist(strsplit(as.character(x %||% ""), ",", fixed = TRUE)))
    entries <- entries[entries != "" & entries != "." & !is.na(entries)]
    if (length(entries) == 0) {
      return(data.frame(gene = character(), region_rank = integer(), stringsAsFactors = FALSE))
    }

    core <- sub(";.*$", "", entries)
    gene <- sub("\\|.*$", "", core)
    region <- ifelse(grepl("\\|", core), sub("^[^|]*\\|", "", core), "")

    region_rank <- ifelse(grepl("exon", region, ignore.case = TRUE), 3L,
                          ifelse(grepl("intron", region, ignore.case = TRUE), 2L,
                                 ifelse(grepl("promoter", region, ignore.case = TRUE), 1L, 0L)))

    out <- data.frame(gene = gene, region_rank = region_rank, stringsAsFactors = FALSE)
    out <- out[out$gene != "" & out$gene != "." & !is.na(out$gene), , drop = FALSE]
    out <- out[order(out$gene, -out$region_rank), , drop = FALSE]
    out[!duplicated(out$gene), , drop = FALSE]
  }

  rec_list <- vector("list", length(keep))
  for (k in seq_along(keep)) {
    i <- keep[[k]]
    g1 <- parse_gene_field(dt$gene1[[i]])
    g2 <- parse_gene_field(dt$gene2[[i]])
    if (nrow(g1) == 0 && nrow(g2) == 0) next

    genes1 <- if (nrow(g1) > 0) g1$gene else character()
    genes2 <- if (nrow(g2) > 0) g2$gene else character()
    both_genes <- intersect(genes1, genes2)

    gg <- bind_rows(g1, g2) %>%
      group_by(gene) %>%
      summarise(region_rank = max(region_rank), .groups = "drop") %>%
      mutate(in_both = gene %in% both_genes)

    rec_list[[k]] <- data.frame(
      sample_id   = rep(sample_id, nrow(gg)),
      gene        = gg$gene,
      event       = rep(event[[i]], nrow(gg)),
      datatype    = "SV",
      eventScore  = evScore[[i]],
      minEv       = minEv[[i]],
      region_rank = gg$region_rank,
      in_both     = gg$in_both,
      stringsAsFactors = FALSE
    )
  }

  rec <- bind_rows(rec_list)
  if (nrow(rec) == 0) return(data.frame())

  rec %>%
    arrange(sample_id, gene,
            desc(eventScore),
            desc(minEv),
            desc(region_rank),
            desc(as.integer(in_both))) %>%
    group_by(sample_id, gene) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    select(sample_id, gene, event, datatype)
}

events_mut_sv <- bind_rows(lapply(seq_len(nrow(manifest)), function(i) {
  sid <- as.character(manifest$sample_id[i])
  bind_rows(
    read_snv_indel_table_events(manifest$snv_path[i], sid),
    read_snv_indel_table_events(manifest$indel_path[i], sid),
    read_sv_events(manifest$sv_path[i], sid)
  )
})) %>%
  filter(!is.na(gene), gene != "", gene != ".",
         !is.na(event), event != "",
         !is.na(datatype), datatype %in% c("MUT", "SV")) %>%
  transmute(
    sample_id = as.character(sample_id),
    gene      = as.character(gene),
    event     = as.character(event),
    datatype  = as.character(datatype)
  ) %>%
  distinct(sample_id, gene, event, datatype)

message("Events (mut+sv): ", nrow(events_mut_sv))

saveRDS(events_mut_sv, file.path(OUT_DIR, "MUT_SV_events.rds"))
fwrite(events_mut_sv, file.path(OUT_DIR, "MUT_SV_events.tsv"), sep = "\t")

### Read snv and indel files
read_snv_indel_variants <- function(path, sample_id) {
  dt <- tryCatch(
    fread(path, sep = "\t", data.table = FALSE, quote = ""),
    error = function(e) {
      warning("Variant read failed: ", basename(path), " | ", conditionMessage(e))
      return(NULL)
    }
  )
  if (is.null(dt) || nrow(dt) == 0) return(data.frame())

  req <- c("#CHROM", "POS", "REF", "ALT", "FILTER", "ANNOVAR_FUNCTION", "EXONIC_CLASSIFICATION")
  miss <- setdiff(req, names(dt))
  if (length(miss) > 0) {
    warning("Variants skip missing columns: ", basename(path), " | ", paste(miss, collapse = ","))
    return(data.frame())
  }

  dt <- dt[is.na(dt$FILTER) | dt$FILTER %in% c("PASS", ".", ""), , drop = FALSE]
  if (nrow(dt) == 0) return(data.frame())

  out <- data.frame(
    sample_id = rep(sample_id, nrow(dt)),
    chrom = norm_chr(as.character(dt[["#CHROM"]])),
    pos   = suppressWarnings(as.integer(dt$POS)),
    ref   = as.character(dt$REF),
    alt   = as.character(dt$ALT),
    func  = as.character(dt$ANNOVAR_FUNCTION),
    exo   = as.character(dt$EXONIC_CLASSIFICATION),
    stringsAsFactors = FALSE
  )

  out <- out[
    !is.na(out$pos) & out$pos > 0 &
      !is.na(out$ref) & !is.na(out$alt) &
      out$chrom %in% std_chr,
    , drop = FALSE
  ]

  distinct(out, sample_id, chrom, pos, ref, alt, .keep_all = TRUE)
}

samples_use <- as.character(manifest$sample_id)

variant_list_all <- lapply(seq_len(nrow(manifest)), function(i) {
  sid <- as.character(manifest$sample_id[i])
  bind_rows(
    read_snv_indel_variants(manifest$snv_path[i], sid),
    read_snv_indel_variants(manifest$indel_path[i], sid)
  )
})

variants_all <- bind_rows(variant_list_all) %>%
  distinct(sample_id, chrom, pos, ref, alt, .keep_all = TRUE)

#### compute TMB
compute_coding_mb_from_gtf <- function(gtf_path) {
  gtf <- rtracklayer::import(gtf_path)
  exons <- gtf[gtf$type == "exon"]
  if (length(exons) == 0) stop("No exons found in GTF: ", gtf_path)
  exons <- reduce(exons, ignore.strand = TRUE)
  sum(width(exons)) / 1e6
}

coding_mb <- compute_coding_mb_from_gtf(gtf_path)

is_nonsyn_for_tmb <- function(func, exo) {
  func <- tolower(func %||% "")
  exo  <- tolower(exo %||% "")
  if (str_detect(exo, "synonymous")) return(FALSE)
  if (str_detect(func, "splic") | str_detect(exo, "splice")) return(TRUE)
  if (str_detect(exo, "frameshift")) return(TRUE)
  if (str_detect(exo, "stopgain|stop_gained|nonsense")) return(TRUE)
  if (str_detect(exo, "nonframeshift|inframe")) return(TRUE)
  if (str_detect(exo, "nonsynonymous|missense")) return(TRUE)
  FALSE
}

compute_tmb_table <- function(variants_df, coding_mb) {
  if (nrow(variants_df) == 0) {
    return(data.frame(sample_id = character(), N_nonsyn = integer(), TMB = numeric(), stringsAsFactors = FALSE))
  }

  keep <- mapply(is_nonsyn_for_tmb, variants_df$func, variants_df$exo, USE.NAMES = FALSE)
  dt <- variants_df[keep, , drop = FALSE]
  tab <- as.data.frame(table(dt$sample_id), stringsAsFactors = FALSE)
  colnames(tab) <- c("sample_id", "N_nonsyn")
  tab$N_nonsyn <- as.integer(tab$N_nonsyn)
  tab$TMB <- tab$N_nonsyn / coding_mb
  tab
}

tmb_df <- compute_tmb_table(variants_all, coding_mb) %>%
  mutate(sample_id = as.character(sample_id))

metadata <- fread(meta_file, sep = "\t", data.table = FALSE)
metadata$sample_id <- as.character(metadata$sample_id)

metadata_with_TMB <- metadata %>%
  left_join(tmb_df %>% select(sample_id, N_nonsyn, TMB), by = "sample_id") %>%
  mutate(
    N_nonsyn = dplyr::coalesce(N_nonsyn, 0L),
    TMB      = dplyr::coalesce(TMB, 0)
  ) %>%
  arrange(factor(sample_id, levels = samples_use))

fwrite(metadata_with_TMB, file.path(OUT_DIR, "metadata_with_TMB.tsv"), sep = "\t")

variants_to_granges_list_snv <- function(variants_df) {
  dt <- variants_df[nchar(variants_df$ref) == 1 & nchar(variants_df$alt) == 1, , drop = FALSE]
  dt <- dt[dt$ref %in% c("A", "C", "G", "T") & dt$alt %in% c("A", "C", "G", "T"), , drop = FALSE]
  if (nrow(dt) == 0) return(NULL)

  dt$chrom_ucsc <- to_ucsc_chr(dt$chrom)
  split_dt <- split(dt, dt$sample_id)

  grl <- lapply(split_dt, function(d) {
    GRanges(
      seqnames = d$chrom_ucsc,
      ranges = IRanges(start = d$pos, end = d$pos),
      REF = d$ref,
      ALT = d$alt
    )
  })
  GRangesList(grl)
}

fit_sbs_signatures <- function(grl, ref_genome) {
  mut_mat <- mut_matrix(vcf_list = grl, ref_genome = ref_genome)

  sigs <- NULL
  for (nm in c("COSMIC_v3.4", "COSMIC_v3.3", "COSMIC_v3.2", "COSMIC_v3")) {
    sigs <- tryCatch(get_known_signatures(nm), error = function(e) NULL)
    if (!is.null(sigs)) break
  }
  if (is.null(sigs)) stop("Could not load known COSMIC signatures via MutationalPatterns.")

  fit_res <- fit_to_signatures(mut_mat, sigs)
  fit_res$contribution
}

if (HAS_MUTPAT && HAS_BSGEN) {
  sig_frac <- NULL
  try({
    grl <- variants_to_granges_list_snv(variants_all)
    if (!is.null(grl)) {
      sig_contrib <- fit_sbs_signatures(grl, BSgenome.Hsapiens.UCSC.hg19)
      sig_frac <- sweep(sig_contrib, 2, colSums(sig_contrib), FUN = "/")
      sig_frac[is.na(sig_frac)] <- 0
    }
  }, silent = TRUE)

  if (!is.null(sig_frac)) {
    sig_df <- as.data.frame(sig_frac)
    sig_df$Signature <- rownames(sig_df)
    sig_df <- sig_df[, c("Signature", setdiff(colnames(sig_df), "Signature")), drop = FALSE]
    fwrite(sig_df, file.path(OUT_DIR, "SBS_signature_fractions.tsv"), sep = "\t")
  } else {
    warning("SBS signature fitting did not run (likely too few SNVs or missing COSMIC signatures).")
  }
}
