# MettleRNASeq: run the blocks in order, or select lines in RStudio.

# Run each stage using the step-by-step manual.

# Replace the temporary output directory with a permanent project folder.

# Stage: load
library(MettleRNASeq)
packageVersion("MettleRNASeq")
expression <- read.csv(system.file("extdata", "mestrado_bone_expression.csv",
  package = "MettleRNASeq"), row.names = 1, check.names = FALSE)
metadata <- read.csv(system.file("extdata", "mestrado_bone_metadata.csv",
  package = "MettleRNASeq"), stringsAsFactors = FALSE)
stopifnot(identical(rownames(expression), metadata$Sample))
dim(expression)
metadata
expression[1:3, 1:4]

# Stage: normalize
set.seed(123)
counts <- as.data.frame(matrix(rnbinom(8 * 60, size = 10, mu = 100), 8))
rownames(counts) <- paste0("Sample_", 1:8)
colnames(counts) <- paste0("Gene_", 1:60)
normalized_synthetic <- normalize_count_matrix(counts, method = "TMM",
  min_count = 10, group_size = 4, remove_low_variance = TRUE)
dim(normalized_synthetic)
attr(normalized_synthetic, "removed_genes")
attr(normalized_synthetic, "removed_samples")
# Continue with the already normalized example matrix; do not normalize it again.
normalized <- expression

# Stage: distance
distances <- compute_distance_matrix(normalized,
  use_correlation = TRUE, correlation_method = "spearman",
  plot_heatmap = TRUE, display_numbers = TRUE,
  main_title = "Sample distance: 1 - Spearman")
round(distances[1:4, 1:4], 3)

# Stage: pca
pca <- pca_analysis(normalized, design_matrix = metadata,
  color_by = "Condition", center = TRUE, scale = FALSE,
  plot_2d = TRUE, plot_3d = FALSE)
data.frame(PC = seq_along(pca$variance),
  variance_percent = pca$variance,
  cumulative_percent = pca$cumulative_variance)
head(pca$cos2)

# Stage: features
pca_selection <- feature_selection_pca(pca$pca_result, normalized,
  num_pcs = 2, num_genes_per_pc = 12,
  select_worst = FALSE, use_last_pcs = FALSE)
length(pca_selection$selected_genes)
head(pca_selection$selected_genes)
refined <- pca_selection$selected_count_matrix
dim(refined)

# Stage: prepare
input <- prepare_input_matrix(refined, design_table = metadata,
  sample_col = "Sample", condition_col = "Condition")
head(input)
attr(input, "sample_map")

# Stage: classification
classification <- classification_approach_advanced(input,
  algorithms = "rf", hyper_params = list(rf = data.frame(mtry = c(1, 3))),
  num_resamples = 3, num_genes = 10, train_percentage = .75,
  metrics = "Accuracy", summary_stat = "median", replace_genes = TRUE,
  parallel = FALSE, seed = 123, include_gene = TRUE,
  error_handling = "stop")
classification$partition_table
head(classification$eval_table)
classification$status_table
head(classification$gene_ids_table)

# Stage: select
scores <- unique(classification$eval_table[,
  c("algorithm", "resample", "summary_test_accuracy")])
scores
chosen_id <- scores$resample[which.max(scores$summary_test_accuracy)]
selection <- genes_from_classification_results(classification$eval_table,
  classification$genes_table, refined, resamples = chosen_id,
  gene_ids_table = classification$gene_ids_table,
  partition_table = classification$partition_table, selection_mode = "unique",
  plot_type = "both", plot_resamples_only = FALSE, return_plot = TRUE)
selection$selected_resamples
selection$selected_gene_ids
selected_expression <- selection$filtered_matrix
dim(selected_expression)

# Stage: acc
stopifnot(identical(metadata$Condition,
  rep(c("FLASH", "STANDARD"), each = 4)))
split <- c(4, 8)
acc <- acc_analysis(selected_expression, condition_split = split,
  method = "spearman", plot = TRUE, orientation = "horizontal")
acc$results

# Stage: icc
icc <- icc_analysis(selected_expression, condition_split = split,
  method = "spearman", plot = TRUE, plot_type = "lower",
  show_coef = TRUE, tl.cex = .6, number.cex = .6, order = "alphabet")
round(icc$correlation_matrix_cond1[1:4, 1:4], 2)
round(icc$correlation_matrix_cond2[1:4, 1:4], 2)

# Stage: differences
differences <- icc_difference(icc, cutoff = NULL,
  unique_pairs = TRUE, plot = TRUE, fontsize = 7)
head(differences$difference_table)
dim(differences$absolute_diff_matrix)

# Stage: difference-filter
large_changes <- icc_difference(icc, cutoff = 1, plot = FALSE)
head(large_changes$difference_table)

# Stage: rules
flash_rules <- mine_gene_associations(selected_expression, condition_rows = 1:4,
  method = "apriori", supp = .5, conf = .8, minlen = 2, maxlen = 2,
  plot_type = "scatter", interactive = FALSE, top_n = 20)
standard_rules <- mine_gene_associations(selected_expression, condition_rows = 5:8,
  method = "apriori", supp = .5, conf = .8, minlen = 2, maxlen = 2,
  plot_type = "scatter", interactive = FALSE, top_n = 20)
length(flash_rules$rules)
length(standard_rules$rules)
head(flash_rules$rules_df)
head(standard_rules$rules_df)
head(arules::itemLabels(flash_rules$transactions))

# Stage: compare
comparison <- compare_gene_rules(flash_rules$rules_df,
  standard_rules$rules_df, condition_names = c("FLASH", "STANDARD"),
  comparison = "genes", layout_type = "linear", plot = TRUE)
comparison$shared_rules
head(comparison$unique_rules_condition1)
head(comparison$unique_rules_condition2)

# Stage: compare-items
item_comparison <- compare_gene_rules(flash_rules$rules_df,
  standard_rules$rules_df, condition_names = c("FLASH", "STANDARD"),
  comparison = "items", plot = FALSE)
nrow(item_comparison$shared_rules)

# Stage: save
directory <- file.path(tempdir(), "MettleRNASeq_step_by_step")
dir.create(directory, recursive = TRUE, showWarnings = FALSE)
saveRDS(list(normalized = normalized, metadata = metadata,
  distances = distances, pca = pca, pca_selection = pca_selection,
  input = input, classification = classification, selection = selection,
  acc = acc, icc = icc, differences = differences,
  flash_rules = flash_rules, standard_rules = standard_rules,
  comparison = comparison), file.path(directory, "full_workflow.rds"))
write.csv(classification$eval_table, file.path(directory, "classification.csv"),
  row.names = FALSE)
write.csv(classification$gene_ids_table, file.path(directory, "sampled_gene_ids.csv"),
  row.names = FALSE)
write.csv(selected_expression, file.path(directory, "selected_expression.csv"))
write.csv(acc$results, file.path(directory, "acc.csv"), row.names = FALSE)
write.csv(icc$correlation_matrix_cond1, file.path(directory, "icc_FLASH.csv"))
write.csv(icc$correlation_matrix_cond2, file.path(directory, "icc_STANDARD.csv"))
write.csv(differences$difference_table, file.path(directory, "icc_differences.csv"),
  row.names = FALSE)
write.csv(flash_rules$rules_df, file.path(directory, "FLASH_rules.csv"),
  row.names = FALSE)
write.csv(standard_rules$rules_df, file.path(directory, "STANDARD_rules.csv"),
  row.names = FALSE)
write.csv(comparison$shared_rules, file.path(directory, "shared_rules.csv"),
  row.names = FALSE)
write.csv(comparison$unique_rules_condition1,
  file.path(directory, "unique_FLASH_rules.csv"), row.names = FALSE)
write.csv(comparison$unique_rules_condition2,
  file.path(directory, "unique_STANDARD_rules.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(directory, "sessionInfo.txt"))
list.files(directory)
