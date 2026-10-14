# Run from R with source(system.file("examples", "run_example.R", package = "MettleRNASeq"))
# Or: Rscript run_example.R output-directory mestrado
run_mettle_example <- function(dataset = c("mestrado", "synthetic"), output_dir = NULL,
                               num_resamples = 3L, plot = TRUE) {
  dataset <- match.arg(dataset)
  if (!requireNamespace("MettleRNASeq", quietly = TRUE)) stop("Install MettleRNASeq first.")
  if (!requireNamespace("randomForest", quietly = TRUE)) stop("Install randomForest for this example.")
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    if (plot) {
      grDevices::pdf(file.path(output_dir, "workflow_plots.pdf"), width = 10, height = 7)
      on.exit(grDevices::dev.off(), add = TRUE)
    }
  }
  set.seed(123)
  if (dataset == "mestrado") {
    expression_file <- system.file("extdata", "mestrado_bone_expression.csv", package = "MettleRNASeq")
    metadata_file <- system.file("extdata", "mestrado_bone_metadata.csv", package = "MettleRNASeq")
    if (!nzchar(expression_file) || !nzchar(metadata_file)) stop("Example data missing; reinstall the documented version.")
    expression <- read.csv(expression_file, row.names = 1, check.names = FALSE)
    metadata <- read.csv(metadata_file, stringsAsFactors = FALSE)
    # These are normalized, preselected MSc values. Do not normalize them again.
    normalized <- expression
  } else {
    metadata <- data.frame(Sample = c(paste0("A", 1:4), paste0("B", 1:4)),
      Condition = rep(c("A", "B"), each = 4), stringsAsFactors = FALSE)
    expression <- as.data.frame(matrix(stats::rnbinom(8 * 60, size = 10, mu = 100), 8))
    rownames(expression) <- metadata$Sample; colnames(expression) <- paste0("Gene_", 1:60)
    expression[metadata$Condition == "A", 1:12] <- expression[metadata$Condition == "A", 1:12] + 60
    normalized <- MettleRNASeq::normalize_count_matrix(expression, method = "TMM", min_count = 10, group_size = 4)
  }
  metadata <- metadata[match(rownames(normalized), metadata$Sample), , drop = FALSE]
  stopifnot(!anyNA(metadata$Sample), identical(rownames(normalized), metadata$Sample))
  distance <- MettleRNASeq::compute_distance_matrix(normalized, plot_heatmap = plot)
  pca <- MettleRNASeq::pca_analysis(normalized, design_matrix = metadata,
    color_by = "Condition", plot_2d = plot)
  features <- MettleRNASeq::feature_selection_pca(pca$pca_result, normalized,
    num_pcs = 2, num_genes_per_pc = 12)
  input <- MettleRNASeq::prepare_input_matrix(features$selected_count_matrix,
    sample_vector = metadata$Sample, condition_vector = metadata$Condition)
  classification <- MettleRNASeq::classification_approach_advanced(input,
    algorithms = "rf", hyper_params = list(rf = data.frame(mtry = c(1, 3))),
    num_resamples = num_resamples, num_genes = 10, seed = 123)
  scores <- unique(classification$eval_table[, c("algorithm", "resample", "summary_test_accuracy")])
  stopifnot(all(is.finite(scores$summary_test_accuracy)))
  best_resample <- scores$resample[which.max(scores$summary_test_accuracy)]
  # Pick one resample to keep the relationship plots small. This is not a significance test.
  selected <- MettleRNASeq::genes_from_classification_results(classification$eval_table,
    classification$genes_table, features$selected_count_matrix, resamples = best_resample,
    partition_table = classification$partition_table, plot_type = "both",
    plot_resamples_only = FALSE, return_plot = plot)
  selected_matrix <- selected$filtered_matrix
  levels <- unique(metadata$Condition)
  n_first <- sum(metadata$Condition == levels[1])
  stopifnot(identical(metadata$Condition, rep(levels, each = n_first)))
  split <- c(n_first, nrow(selected_matrix))
  acc <- MettleRNASeq::acc_analysis(selected_matrix, condition_split = split, plot = plot)
  acc_table <- if (plot) acc$results else acc
  icc <- MettleRNASeq::icc_analysis(selected_matrix, condition_split = split, plot = plot,
    order = "alphabet", tl.cex = .6)
  difference <- MettleRNASeq::icc_difference(icc, plot = plot, fontsize = 6)
  # Numeric columns use arules' discretization. See the manual for explicit state definitions.
  first <- MettleRNASeq::mine_gene_associations(selected_matrix, seq_len(n_first),
    supp = .5, conf = .8, maxlen = 2, interactive = FALSE, top_n = 20)
  second <- MettleRNASeq::mine_gene_associations(selected_matrix, seq.int(n_first + 1L, nrow(selected_matrix)),
    supp = .5, conf = .8, maxlen = 2, interactive = FALSE, top_n = 20)
  rules <- MettleRNASeq::compare_gene_rules(first$rules_df, second$rules_df,
    condition_names = levels, comparison = "genes", plot = plot)
  summary <- data.frame(dataset = dataset, input_samples = nrow(expression),
    input_genes = ncol(expression), pca_genes = length(features$selected_genes),
    resamples = num_resamples, partitions = nrow(classification$partition_table),
    selected_resample = best_resample, selected_genes = length(selected$selected_gene_ids),
    rules_condition1 = length(first$rules), rules_condition2 = length(second$rules),
    shared_gene_rules = nrow(rules$shared_rules))
  if (!is.null(output_dir)) {
    export <- function(x, name) utils::write.csv(x, file.path(output_dir, paste0(name, ".csv")), row.names = FALSE)
    export(summary, "summary"); export(metadata, "metadata")
    export(classification$eval_table, "classification"); export(classification$partition_table, "partitions")
    export(classification$sample_map, "replicate_map"); export(classification$status_table, "fit_status")
    export(classification$gene_ids_table, "resampled_gene_ids")
    export(classification$test_eval_table, "heldout_replicate_accuracy")
    export(data.frame(Gene = selected$selected_gene_ids), "selected_genes")
    export(acc_table, "across_condition_correlations"); export(difference$difference_table, "correlation_differences")
    export(first$rules_df, "rules_condition1"); export(second$rules_df, "rules_condition2")
    export(rules$shared_rules, "shared_gene_rules")
    utils::capture.output(utils::sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
    saveRDS(list(summary = summary, classification = classification$eval_table,
      gene_ids = classification$gene_ids_table, sample_map = classification$sample_map,
      selected_gene_ids = selected$selected_gene_ids, selected_matrix = selected_matrix,
      icc = icc, rules = list(first = first$rules, second = second$rules)),
      file.path(output_dir, "example_results.rds"))
  }
  print(summary)
  invisible(list(summary = summary, classification = classification,
    selected = selected, pca = pca, distance = distance, rules = rules))
}

if (sys.nframe() == 0L) {
  arguments <- commandArgs(trailingOnly = TRUE)
  directory <- if (length(arguments)) arguments[1] else "MettleRNASeq_example"
  dataset <- if (length(arguments) >= 2L) arguments[2] else "mestrado"
  run_mettle_example(dataset, output_dir = directory)
}
