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
  accession numbers and sources, are listed in Supplementary Table X of the
  manuscript.
- **Proteomics data** generated in this study are available from
  ProteomeXchange/PRIDE under accession **PXD[xxxxxx]**.
- **In-house RNA-seq and WGS data** are available [under controlled access from
  EGA, accession EGAS00001004813].


To rerun the analyses, place the data in the folders below (file names must
match those used in the scripts).

### Discovery cohort (RNA-seq) – raw data folder

```
<RNA-seq raw data folder>/
├── meta/raw_metadata/aspera/   # in-house sample sheet + clinical information (S0)
├── meta/raw_metadata/NCBI/     # SRA run tables (*.csv) of the public BioProjects (S1)
├── data/raw_featurecount/<dataset>/   # featureCounts output (*.biotype_counts_mqc.tsv) per dataset (S2)
└── data/rds/                   # count matrices obtained from publications/GEO (S4)
```

TCGA-SARC data are downloaded by `S3_TCGA_dataset/S1_download_raw_tcga_data.R`.

### Validation cohort (microarray) – raw data folder

All five datasets are public (GEO: GSE21050, GSE21122, GSE30929, GSE71118,
GSE159659). The folders below are created by the scripts, except
`meta/raw_metadata/`:

```
<microarray raw data folder>/
├── meta/raw_metadata/          # GEO series matrix files, GSExxxxx_series_matrix.txt (download manually, unzipped) (S2)
├── meta/lps_metadata/          # liposarcoma / normal fat samples per dataset (S3)
├── meta/metadata/              # integrated metadata (S8, S11)
├── raw_data/<GSE>/<GSE>_RAW/   # CEL files downloaded from GEO (S1)
├── data/lps_data/<GSE>_LPS/    # CEL files of selected samples (S4)
├── qc_reports/                 # arrayQualityMetrics reports (S5)
└── data/rds/                   # gene-level expression per dataset and integrated matrices (S6–S11)
```

### WGS – raw data folder

```
<WGS raw data folder>/
├── WGS_input/snv/somatic_functional/    # somatic SNV calls per sample
├── WGS_input/indels/somatic_functional/ # somatic indel calls per sample
├── WGS_input/sophia/high_conf/          # structural variants (SOPHIA)
├── WGS_input/aceseq_cnv/                # copy-number calls (ACEseq)
├── WGS_input/gencode.v19.annotation.gtf # gene annotation (GENCODE v19, hg19)
├── WGS_input/<sample sheet>.tsv, <clinical information>.xlsx   # in-house metadata (S2)
└── WGS_input/renamed/, metadata/, preprocessed/   # created by S1–S5
```

Results are written to `WGS_results/` (figures in `WGS_results/Figures/`).

### Analysis modules – `Inputs/` folders

| Module | Required files in `Inputs/` | Origin |
|---|---|---|
| Discovery_cohort_RNAseq | `Batch_corrected_data_discovery_all.rds`, `Batch_corrected_data_discovery_WD_DD.rds`, `Metadata_discovery_all.csv`, `Metadata_discovery_WD_DD.csv`, `aspera_rna_metadata_with_survival.csv`, `TCGA_sarc_dd_survival_data.csv` | Output of `01_data_integration/` S0–S6 (see table below) |
| | `GSEA_custom_signatures/TLS_Sautes-Fridman_2019.csv` | TLS gene signature from Sautès-Fridman *et al.*, *Nat Rev Cancer* 2019 |
| | `DEA_validation_files/Cluster_[1-6]_DEGs.csv` | Copy of `Validation_cohort_microarray/Outputs/Files/Limma/Cluster_[1-6]_DEGs_all.csv` |
| Validation_cohort_microarray | `Batch_corrected_data_validation_all.rds`, `Batch_corrected_data_validation_WD_DD.rds`, `Metadata_validation_all.csv`, `Metadata_validation_WD_DD.csv` | Output of `01_data_integration/` S1–S11 (see table below) |
| | `survival_data_GSE21050.csv`, `survival_data_GSE71118.csv` | Metastasis-free survival of GSE21050 / GSE71118 [source to be added] |
| Proteomics | `Project_968_878_…_TJ.xlsx`, `Project_1055_…_TJ.xlsx`, `Project_1172_…_TJ.xlsx` (DIA-NN protein reports), `Proteomic_metadata.xlsx` | Proteomics data of this study |
| | `Transcriptomic_files/Cluster_[1-6]_DEGs_all.csv` | Copy of `Discovery_cohort_RNAseq/Outputs/Files/Limma/Cluster_[1-6]_DEGs_all.csv` |

The final outputs of the data-integration scripts have the same file names as
the files read by `02_analysis_and_figure_scripts/`; copy them into the
module's `Inputs/` folder.

Discovery cohort (`Discovery_cohort_RNAseq/Inputs/`):

| File | Written by |
|---|---|
| `Batch_corrected_data_discovery_all.rds` | S5 (`data/rds/`) |
| `Batch_corrected_data_discovery_WD_DD.rds` | S6 (`data/WD_DD/`) |
| `Metadata_discovery_all.csv` | S4 (`meta/metadata/`) |
| `Metadata_discovery_WD_DD.csv` | S6 (`meta/WD_DD/`) |
| `aspera_rna_metadata_with_survival.csv` | S0 (`meta/metadata/`) |
| `TCGA_sarc_dd_survival_data.csv` | `S3_TCGA_dataset/S2` (`meta/metadata/`) |

Validation cohort (`Validation_cohort_microarray/Inputs/`):

| File | Written by |
|---|---|
| `Batch_corrected_data_validation_all.rds` | S10 (`data/rds/`) |
| `Batch_corrected_data_validation_WD_DD.rds` | S11 (`data/rds/`) |
| `Metadata_validation_all.csv` | S8 (`meta/metadata/`) |
| `Metadata_validation_WD_DD.csv` | S11 (`meta/metadata/`) |

---

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

### Script overview and figure panels

Scripts listed without a folder prefix are in the module's
`02_analysis_and_figure_scripts/` folder. Figure panels refer to the published
article. Panels 1G, 3G, 5H and Figure 6 are schematic diagrams and panel S1H
shows H&E images; these were not generated with code.

<details>
<summary><b>Transcriptomics – Discovery cohort (RNA-seq)</b></summary>

| Script | Description | Figure |
|---|---|---|
| `01_data_integration/S0_aspera_metadata.R` | In-house RNA-seq sample and clinical metadata | – |
| `01_data_integration/S1_fastq_metadata.R` | Select liposarcoma / normal fat samples from SRA run tables | – |
| `01_data_integration/S2_featurecount_combination.R` | Combine featureCounts output; Ensembl → gene symbol (biomaRt) | – |
| `01_data_integration/S3_TCGA_dataset/S1–S3` | Download TCGA-SARC, select liposarcomas, extract survival, gene symbols | – |
| `01_data_integration/S4_integrate_rna_count_data.R` | Merge count matrices and metadata of all datasets | – |
| `01_data_integration/S5_normalization_batch_correction.R` | TMM + voom normalisation, ComBat batch correction (all samples) | – |
| `01_data_integration/S6_WD_DD_normalization_batch_correction.R` | Same, WD/DD liposarcoma subset (n = 335) | – |
| `1_Treemap_1A.R` | Distribution of datasets | 1A |
| `2_Barplot_datasets_S1A.R` | Histological subtypes per dataset | S1A |
| `3_UMAP_all_samples_1C_S1C_2A.R` | UMAP of all samples by subtype and by dataset; WD/DD subset highlighted | 1C, S1C, 2A (top) |
| `4_UMAP_expression_1D.R` | UMAP coloured by *MDM2*, *PPARG* and *TOP2A* expression | 1D |
| `5_Consensus_clustering_S2AB.R` | Consensus clustering (k-means, 1,000 resamplings); delta area plot and consensus matrix | S2A, S2B |
| `6_UMAP_clusters_2A.R` | UMAP of the WD/DD subset by cluster | 2A (bottom) |
| `7_DEA_clusters.R` | limma: each cluster vs. all others | – |
| `8_Cluster_heatmap_2B.R` | Heatmap of cluster marker genes | 2B |
| `9_Cluster_correlation_RNA_vs_RNA_S2E.R` | Correlation between clusters within the discovery cohort | S2E |
| `10_Cluster_correlation_RNA_vs_MA_2D_S3.R` | Correlation between discovery and validation clusters | 2D, S3 |
| `11_ssGSEA.R` | ssGSEA (GSVA) scores of selected gene sets | – |
| `12_ssGSEA_heatmap_2E.R` | Heatmap of mean ssGSEA scores per cluster | 2E |
| `13_Survival_RNA_2F.R` | Overall survival by cluster | 2F |
| `14_Immune_infiltration_violins_3D_S5C.R` | Immune cell fractions (EPIC) and pairwise tests (Dunn) | 3D, S5C |
| `15_DEA_clusters_subsets.R` | limma: cluster 1 vs. cluster 2 | – |
| `16_Immune_marker_heatmap_3E.R` | Immune marker expression in clusters 1 and 2 | 3E |
| `17_GSEA_cluster1_S5D.R` | GSEA in cluster 1 (Hallmark IFN-γ response, TLS signature) | S5D |
| `18_DEA_volcano_clusters.R` | Volcano plots | 4A, 5A |
| `19_ORA_clusters_RNA.R` | Over-representation analysis per cluster (tables) | – |
| `20_ORA_replotting_4B_5B.R` | Over-representation analysis bar plots | 4B, 5B |
| `21_Dotplot_4C.R` | TGF-β / MSC-related gene expression per cluster | 4C |
| `22_Survival_gene_expression_4EF_5G.R` | Overall survival by *CD34*, *TGFBR3* and *PROM1* expression (non-adipocytic tumours) | 4E, 4F, 5G |
| `23_GSEA_cluster4_S6A.R` | GSEA in cluster 4 (Wnt and stem cell gene sets) | S6A |
| `24_TF_analysis_violin_5C.R` | GLI1 transcription factor activity (decoupleR/VIPER, CollecTRI) | 5C |
| `25_Violin_plots_clusters_5F_S6G.R` | Expression of *ACAN*, *COL2A1*, *PROM1*, *SOX2*, *SOX9* and *RUNX2* per cluster | 5F, S6G |
| `pathway_analysis_function.R` | Helper function used by script 19 | – |

</details>

<details>
<summary><b>Transcriptomics – Validation cohort (microarray)</b></summary>

| Script | Description | Figure |
|---|---|---|
| `01_data_integration/S1_download_microarray_CEL_file.R` | Download raw CEL files from GEO | – |
| `01_data_integration/S2_handle_gse_metadata.R` | Parse sample annotation from GEO series matrix files | – |
| `01_data_integration/S3_extract_lps_samples.R` | Select liposarcoma and normal fat samples | – |
| `01_data_integration/S4_select_lps_raw_data.R` | Collect CEL files of the selected samples | – |
| `01_data_integration/S5_raw_microarray_QC.R` | Array QC before and after RMA (arrayQualityMetrics) | – |
| `01_data_integration/S6_GSE159659.R` | GSE159659 (Clariom S): normalised data from GEO, gene symbols | – |
| `01_data_integration/S7_microarray_to_gene_expr_bioconductor.R` | RMA normalisation; probe → gene symbol (Affymetrix HG-U133A / Plus 2.0) | – |
| `01_data_integration/S8_combine_microarray_metadata.R` | Combine metadata, harmonise subtypes, remove arrays failing QC | – |
| `01_data_integration/S9_combine_microarray_exprs_data.R` | Collapse probes to genes (mean), merge datasets on shared genes, remove arrays failing QC | – |
| `01_data_integration/S10_microarray_remove_batch_effect.R` | ComBat batch correction (all samples) | – |
| `01_data_integration/S11_WD_DD_normalization_batch_correction.R` | Same, WDLPS/DDLPS subset (n = 270) | – |
| `1_Treemap_1B.R` | Distribution of datasets | 1B |
| `2_Barplot_datasets_S1B.R` | Histological subtypes per dataset | S1B |
| `3_UMAP_all_samples_1E_S1D.R` | UMAP of all samples by subtype and by dataset | 1E, S1D |
| `4_UMAP_expression_S1E.R` | UMAP coloured by *MDM2*, *PPARG* and *TOP2A* expression | S1E |
| `5_Consensus_clustering_S2CD.R` | Consensus clustering; delta area plot and consensus matrix | S2C, S2D |
| `6_DEA_clusters_MA.R` | limma: each cluster vs. all others | – |
| `7_Cluster_heatmap_2C.R` | Heatmap of cluster marker genes | 2C |
| `8_Cluster_correlation_MA_vs_MA_S2F.R` | Correlation between clusters within the validation cohort | S2F |
| `9_ssGSEA.R` | ssGSEA (GSVA) scores of selected gene sets | – |
| `10_ssGSEA_heatmap_S2G.R` | Heatmap of mean ssGSEA scores per cluster | S2G |
| `11_Survival_metastasis_MA_2G.R` | Metastasis-free survival by cluster | 2G |

</details>

<details>
<summary><b>Proteomics</b></summary>

| Script | Description | Figure |
|---|---|---|
| `01_data_integration/1_Proteomic_data_integration.R` | Merge DIA-NN batches; QC, normalisation, imputation | – |
| `01_data_integration/2_Batch_correction.R` | ComBat batch correction | – |
| `01_data_integration/3_Raw_data_merge.R` | Merge non-imputed data (protein detection analyses) | – |
| `4_UMAP_all_1H_S1F_S1G.R` | UMAP of all samples; CDK4, ADIPOQ and TOP2B expression; multi-sample patients | 1H, S1F, S1G |
| `5_Consensus_clustering_tumor_S5A.R` | Consensus clustering of tumour samples | S5A |
| `6_UMAP_clusters_detection_3A_S6CF.R` | UMAP by cluster; GLI1 and HHIP detection (UMAP, bar plots, Fisher's exact test) | 3A, S6C–F |
| `7_Barplot_cluster_distribution_S5B.R` | Cluster distribution of samples from the same tumour | S5B |
| `8_Limma_clusters_tumor.R` | limma: each cluster vs. all others | – |
| `9_Correlation_RNA_protein_3B_S4.R` | Correlation between discovery (RNA) and proteomic clusters | 3B, S4 |
| `10_ssGSEA_proteomics.R` | ssGSEA (GSVA) scores of selected gene sets | – |
| `11_ssGSEA_heatmap_proteomics_3C.R` | Heatmap of mean ssGSEA scores per cluster | 3C |
| `12_Dotplot_cluster_markers_3F_4D.R` | Immune and TGF-β / MSC-related protein expression per cluster | 3F, 4D |

</details>

<details>
<summary><b>WGS</b></summary>

| Script | Description | Figure |
|---|---|---|
| `01_data_processing/S1_handle_file_names.R` | Harmonise sample names; keep samples with SNV, indel, SV and CNV calls | – |
| `01_data_processing/S2_dna_metadata.R` | DNA sample and clinical metadata | – |
| `01_data_processing/S3_oncoprint_manifest.R` | File manifest per sample | – |
| `01_data_processing/S4_oncoprint_snv_sv_processing.R` | SNV/indel/SV events, TMB, SBS signatures | – |
| `01_data_processing/S5_oncoprint_CNA_processing.R` | Gene-level copy-number events (ACEseq) | – |
| `01_data_processing/S6_oncoprint_plot_1F.R` | OncoPrint of samples with matched genomic and transcriptomic data (hg19) | 1F |
| `02_analysis_and_figure_scripts/1_CNV_oncoprint_GLI1_5D_5E.R` | Copy-number OncoPrint ordered by cluster; *GLI1* amplification by cluster (Fisher's exact test) | 5D, 5E |

`02_analysis_and_figure_scripts/1_CNV_oncoprint_GLI1_5D_5E.R` additionally needs, in the raw WGS data folder:
`WGS_input/preprocessed/data_cna.txt` (TCGA-SARC discrete copy-number calls, downloaded from cBioPortal [study and version to be added]),
`WGS_input/preprocessed/aspera_final_rna_metadata.csv` (output of RNA-seq S0) and
`clusters_k6.csv` (copy of `Transcriptomics/Discovery_cohort_RNAseq/Outputs/Files/Consensus_clustering/Consensus_clustering.csv`).



---

