# Optional larger search on the packaged processed-expression example.
# Install the development version containing summarize_gene_resampling().
library(MettleRNASeq)
expression <- read.csv(system.file("extdata", "bone_expression.csv",
  package = "MettleRNASeq"), row.names = 1, check.names = FALSE)
metadata <- read.csv(system.file("extdata", "bone_metadata.csv",
  package = "MettleRNASeq"), stringsAsFactors = FALSE)
input <- prepare_input_matrix(expression, sample_vector = metadata$Sample,
  condition_vector = metadata$Condition)
classification <- classification_approach_advanced(input, algorithms = "rf",
  hyper_params = NULL, tune_length = 10, num_resamples = 100,
  num_genes = 10, seed = 123, store_models = FALSE)
diagnostics <- summarize_gene_resampling(classification, cutoff = .85)
print(diagnostics$coverage_table)
# Select the union of genes from all qualifying sets. Empty selection is valid.
selection <- NULL
if (diagnostics$coverage_table$qualifying_resamples > 0) {
  selection <- genes_from_classification_results(classification$eval_table,
    classification$genes_table, expression, cutoff = .85,
    metric_column = "summary_test_accuracy", algorithm = "rf",
    selection_mode = "unique", return_plot = FALSE)
  print(selection$selected_gene_ids)
} else {
  message("No gene sets reached the specified cutoff; inspect coverage and scores.")
}
saveRDS(list(classification = classification, diagnostics = diagnostics,
  selection = selection, session = sessionInfo()), "resampling_diagnostics.rds")
