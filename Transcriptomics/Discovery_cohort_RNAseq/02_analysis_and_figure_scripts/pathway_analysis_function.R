# Over-representation analysis with clusterProfiler::enricher().
# Gene sets with BH-adjusted p < 0.05 are written to csv_file; if no gene set passes,
# all results are written to "<csv_file>_allResults.csv" instead.
pathway_analysis <- function(gene_list,
                             gene_sets,
                             background_genes,
                             csv_file) {
  enrich_res <- enricher(gene = gene_list,
                         universe = background_genes,
                         TERM2GENE = gene_sets,
                         pvalueCutoff = 0.05,
                         pAdjustMethod = "BH",
                         minGSSize = 5,
                         maxGSSize = 500,
                         qvalueCutoff = 0.05)

  if (is.null(enrich_res) || nrow(enrich_res@result) == 0) {
    message("No enriched terms found for: ", csv_file)
    return(invisible(NULL))
  }

  enrich_table <- enrich_res@result
  enrich_table_filtered <- enrich_table[enrich_table$p.adjust < 0.05, ]

  if (nrow(enrich_table_filtered) > 0) {
    write_csv(enrich_table_filtered, file = csv_file)
  } else {
    message("No terms passed the p.adjust < 0.05 filter for: ", csv_file)
    write_csv(enrich_table, file = sub(".csv", "_allResults.csv", csv_file))
  }
  invisible(enrich_res)
}
