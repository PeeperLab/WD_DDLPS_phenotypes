# Molecular subtypes of well-differentiated / dedifferentiated liposarcoma

Analysis code accompanying the manuscript:

This repository contains the R code used to integrate and analyse the bulk
transcriptomic (RNA-seq discovery cohort, microarray validation cohort),
proteomic and whole-genome sequencing (WGS) data in the study, and to generate
the figures. **It contains code only; no data are included** (see
[Data availability].

---

## Repository structure

```
WD_DDLPS_phenotypes/
├── Transcriptomics/
│   ├── Discovery_cohort_RNAseq/
│   │   ├── 01_data_integration/             # S0–S6: build the integrated, batch-corrected RNA-seq cohort
│   │   │   └── S3_TCGA_dataset/             #   download and process TCGA-SARC cohort
│   │   └── 02_analysis_and_figure_scripts/  # 1–25: analyses and figures
│   └── Validation_cohort_microarray/
│       ├── 01_data_integration/             # S1–S11: build the integrated, batch-corrected microarray cohort
│       └── 02_analysis_and_figure_scripts/  # 1–11: analyses and figures
├── Proteomics/
│   ├── 01_data_integration/                 # 1–3: merge batches, QC, imputation, batch correction
│   └── 02_analysis_and_figure_scripts/      # 4–12: analyses and figures
├── WGS/
│   ├── 01_data_processing/                  # S1–S5: in-house WGS processing; S6: OncoPrint (Fig. 1F)
│   └── 02_analysis_and_figure_scripts/      # CNV OncoPrint by cluster, GLI1 amplification (Fig. 5D, 5E)
├── renv.lock                                # exact R package versions (see "Software requirements")
├── .Rprofile, renv/                         # renv set-up (activates the project library)
├── .gitignore                               # whitelist: only code/docs can be committed
└── README.md
```

In each module, `01_*` prepares the data and `02_analysis_and_figure_scripts`
contains the code for all analyses and figure panels in the manuscript.

`Proteomics/`, `Transcriptomics/Discovery_cohort_RNAseq/` and
`Transcriptomics/Validation_cohort_microarray/` use the same layout at run time:

```
<module>/
├── 01_data_integration/             # scripts (in this repository)
├── 02_analysis_and_figure_scripts/  # scripts (in this repository)
├── Inputs/                          # input data – NOT in this repository
└── Outputs/                         # created automatically by the scripts
    ├── Figures/
    └── Files/
```

---

## Software requirements

- R 4.5.2
- Exact versions of all R packages are recorded in `renv.lock`. To install
  them, open the repository folder in R (or the `.Rproj` file in RStudio) and run:

  ```r
  install.packages("renv")
  renv::restore()
  ```

  Alternatively, install the packages listed below manually.

Main packages:

| Purpose | Packages |
|---|---|
| Data handling | tidyverse (dplyr, tidyr, readr, stringr, tibble, ggplot2), data.table, readxl, openxlsx, reshape2, matrixStats |
| Normalisation / batch correction | edgeR, limma, sva (ComBat) |
| Microarray processing | affy (RMA), arrayQualityMetrics, GEOquery, hgu133plus2.db, hgu133a.db |
| Dimensionality reduction / clustering | umap, ConsensusClusterPlus |
| Differential expression / enrichment | limma, GSVA, msigdbr, clusterProfiler, fgsea, enrichplot |
| Deconvolution / TF activity | EPIC, decoupleR, OmnipathR, viper |
| Survival | survival, survminer |
| Visualisation | ComplexHeatmap, circlize, pheatmap, ggpubr, ggrepel, patchwork, treemapify, RColorBrewer, scales, ggraph, tidygraph |
| Genomics (WGS) | GenomicRanges, IRanges, rtracklayer, MutationalPatterns, BSgenome.Hsapiens.UCSC.hg19 |
| Public data access | TCGAbiolinks, SummarizedExperiment, biomaRt |

Installation:

```r
install.packages(c("BiocManager", "tidyverse", "data.table", "readxl", "openxlsx",
                   "reshape2", "matrixStats", "umap", "msigdbr", "survival", "survminer",
                   "pheatmap", "ggpubr", "ggrepel", "patchwork", "treemapify",
                   "RColorBrewer", "scales", "ggraph", "tidygraph", "rstatix", "circlize"))
BiocManager::install(c("edgeR", "limma", "sva", "ConsensusClusterPlus", "GSVA",
                       "clusterProfiler", "fgsea", "enrichplot", "decoupleR", "OmnipathR",
                       "viper", "ComplexHeatmap", "GenomicRanges", "IRanges", "rtracklayer",
                       "MutationalPatterns", "BSgenome.Hsapiens.UCSC.hg19",
                       "TCGAbiolinks", "SummarizedExperiment", "biomaRt",
                       "affy", "arrayQualityMetrics", "GEOquery", "hgu133plus2.db", "hgu133a.db"))
remotes::install_github("GfellerLab/EPIC")
```

Some scripts need internet access at run time: `biomaRt` (Ensembl gene
annotation), `msigdbr` (MSigDB gene sets), `OmnipathR` (CollecTRI regulons) and
`TCGAbiolinks` (TCGA-SARC download) and `GEOquery` (GEO microarray data). Results can differ slightly with newer
versions of these online resources.

---

## Setting the working directory

All file paths inside the scripts are relative, so **before running
a script, set the working directory** (e.g. with `setwd()` or the RStudio
*Session → Set Working Directory* menu) to the folder given in the comment at
the top of that script:

| Scripts | Working directory |
|---|---|
| `Proteomics/01_data_integration/*`, `Proteomics/02_analysis_and_figure_scripts/*` | `<repository>/Proteomics` |
| `Transcriptomics/Discovery_cohort_RNAseq/02_analysis_and_figure_scripts/*` | `<repository>/Transcriptomics/Discovery_cohort_RNAseq` |
| `Transcriptomics/Validation_cohort_microarray/02_analysis_and_figure_scripts/*` | `<repository>/Transcriptomics/Validation_cohort_microarray` |
| `Transcriptomics/Discovery_cohort_RNAseq/01_data_integration/*` | your folder with the raw RNA-seq data (layout below) |
| `Transcriptomics/Validation_cohort_microarray/01_data_integration/*` | your folder with the raw microarray data (layout below) |
| `WGS/01_data_processing/*`, `WGS/02_analysis_and_figure_scripts/*` | your folder with the raw WGS data (layout below) |

Example:

```r
setwd("/path/to/WD_DDLPS_phenotypes/Proteomics")
source("02_analysis_and_figure_scripts/4_UMAP_all_1H_S1F_S1G.R")
```

---

## Data availability

No data are included in this repository.

- **Public datasets** used in this study (RNA-seq, microarray, TCGA-SARC), with
  accession numbers and sources, are listed in Supplementary Table of the
  manuscript.
- **Proteomics data** generated in this study are available from
  ProteomeXchange/PRIDE under accession **PXD[xxxxxx]**.
- **In-house RNA-seq and WGS data** are available [under controlled access from
  EGA, accession EGASEGAS00001004813].


## Running the analyses

Run the scripts in numerical order within each folder. Because some modules use
results from other modules, the overall order is:

1. **Discovery cohort – data integration**: `Discovery_cohort_RNAseq/01_data_integration/` S0 → S1 → S2 → `S3_TCGA_dataset/` S1–S3 → S4 → S5 → S6; copy the outputs into `Discovery_cohort_RNAseq/Inputs/` (same file names, table above).
2. **Discovery cohort – analyses**: `Discovery_cohort_RNAseq/02_analysis_and_figure_scripts/` 1–9.
3. **Validation cohort**: `Validation_cohort_microarray/01_data_integration/` S1 → S11; copy the outputs into `Validation_cohort_microarray/Inputs/` (same file names, table above); run `02_analysis_and_figure_scripts/` 1–11; copy the cluster DEG tables into `Discovery_cohort_RNAseq/Inputs/DEA_validation_files/`.
4. **Discovery cohort – analyses**: `Discovery_cohort_RNAseq/02_analysis_and_figure_scripts/` 10–25.
5. **Proteomics**: `Proteomics/01_data_integration/` 1–3; copy the discovery DEG tables into `Proteomics/Inputs/Transcriptomic_files/`; run `Proteomics/02_analysis_and_figure_scripts/` 4–12.
6. **WGS**: `WGS/01_data_processing/` S1–S6, then `WGS/02_analysis_and_figure_scripts/1_CNV_oncoprint_GLI1_5D_5E.R`.

Random seeds are fixed in the scripts (`set.seed(1507)` for clustering, UMAP and
enrichment analyses; `set.seed(123)` for data integration and WGS).


