# Complete real-data workflow: keep the 235 previously selected candidate genes.
# Automatic RF tuning; union of all maximum-score sets; no fixed score cutoff.
source(system.file("examples", "run_example.R", package = "MettleRNASeq"))
result <- run_mettle_example("bone", output_dir = "bone_results",
  search_mode = "expanded", num_resamples = 100, plot = FALSE)
result$diagnostics$coverage_table
unique(result$classification$eval_table[, c("resample", "summary_test_accuracy")])
result$selected$selected_gene_ids

# Choose num_resamples explicitly (e.g. 10000) for a more extensive search.
# More resampling increases runtime and may change the selected gene sets.
# Relative best-score selection is exploratory and does not establish significance.
