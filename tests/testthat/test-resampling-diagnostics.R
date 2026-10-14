diagnostic_input <- function() {
  set.seed(42)
  x <- matrix(rpois(8 * 8, 30), 8, 8)
  x[5:8, ] <- x[5:8, ] + 20
  dimnames(x) <- list(paste0("sample", 1:8), paste0("gene", 1:8))
  prepare_input_matrix(as.data.frame(x), sample_vector = rownames(x),
    condition_vector = rep(c("A", "B"), each = 4))
}

test_that("automatic tuning and explicit equivalent grids have the same scores", {
  input <- diagnostic_input()
  automatic <- classification_approach_advanced(input, "rf", num_resamples = 1,
    num_genes = 4, tune_length = 2)
  grid <- automatic$models$rf$resample_1_split_1$results[, "mtry", drop = FALSE]
  explicit <- classification_approach_advanced(input, "rf", list(rf = grid),
    num_resamples = 1, num_genes = 4)
  expect_equal(automatic$eval_table, explicit$eval_table)
  expect_equal(automatic$gene_ids_table, explicit$gene_ids_table)
  expect_equal(automatic$tuning_table, explicit$tuning_table)
  expect_identical(automatic$configuration$tuning_strategy, "automatic")
  expect_identical(explicit$configuration$tuning_strategy, "explicit_grid")
})

test_that("discarding models preserves all scores, genes and tuning decisions", {
  input <- diagnostic_input()
  stored <- classification_approach_advanced(input, "rf", list(rf = data.frame(mtry = 1)), 2, 4)
  light <- classification_approach_advanced(input, "rf", list(rf = data.frame(mtry = 1)), 2, 4,
    store_models = FALSE)
  for (name in c("eval_table", "partition_table", "genes_table", "gene_ids_table",
                 "test_eval_table", "status_table", "tuning_table")) {
    expect_equal(stored[[name]], light[[name]])
  }
  expect_length(light$models$rf, 0L)
  expect_length(stored$models$rf, 8L)
  expect_lt(as.numeric(object.size(light)), as.numeric(object.size(stored)))
  expect_error(classification_approach_advanced(input, "rf", num_resamples = 1,
    num_genes = 4, tune_length = 0), "tune_length")
  expect_error(classification_approach_advanced(input, "rf", num_resamples = 1,
    num_genes = 4, store_models = NA), "Logical flags")
})

resampling_fixture <- function() {
  draws <- data.frame(Resample_1 = c("g1", "g1", "g2"),
    Resample_2 = c("g2", "g3", "g3"), Resample_3 = c("g1", "g3", "g3"))
  legacy <- data.frame(Resample_1 = c(1, 1, 2), Resample_2 = c(2, 3, 3), Resample_3 = c(1, 3, 3))
  attr(legacy, "gene_order") <- c("g1", "g2", "g3", "g4")
  scores <- data.frame(algorithm = rep(c("rf", "ranger"), each = 6),
    resample = rep(rep(1:3, each = 2), 2),
    summary_test_accuracy = rep(c(.9, .5, NA, .8, .95, .6), each = 2),
    status = rep(c("ok", "ok", "error", "ok", "ok", "ok"), each = 2))
  list(eval_table = scores, genes_table = legacy, gene_ids_table = draws)
}

test_that("coverage includes unseen genes and recurrence excludes duplicate draws", {
  result <- summarize_gene_resampling(resampling_fixture(), cutoff = .8)
  rf <- result$gene_table[result$gene_table$algorithm == "rf", ]
  expect_identical(rf$gene_id, c("g1", "g2", "g3", "g4"))
  expect_identical(rf$resamples_sampled, c(2L, 2L, 2L, 0L))
  expect_identical(rf$resamples_evaluated, c(1L, 2L, 1L, 0L))
  expect_identical(rf$resamples_qualified, c(1L, 1L, 0L, 0L))
  expect_equal(rf$qualification_fraction, c(1, .5, 0, NA))
  expect_equal(rf$qualified_presence_fraction, c(1, 1, 0, 0))
  coverage <- result$coverage_table[result$coverage_table$algorithm == "rf", ]
  expect_equal(coverage$sampling_coverage, .75)
  expect_equal(coverage$evaluated_resamples, 2)
  expect_equal(coverage$qualifying_resamples, 1)
  expect_equal(coverage$min_resamples_per_gene, 0)
  ranger <- result$gene_table[result$gene_table$algorithm == "ranger", ]
  expect_identical(ranger$resamples_qualified, c(1L, 2L, 1L, 0L))
})

test_that("undefined qualification and empty qualifying sets remain explicit", {
  no_cutoff <- summarize_gene_resampling(resampling_fixture(), algorithm = "rf")
  expect_true(all(is.na(no_cutoff$gene_table$resamples_qualified)))
  expect_true(all(is.na(no_cutoff$gene_table$qualification_fraction)))
  none <- summarize_gene_resampling(resampling_fixture(), cutoff = 1, algorithm = "rf")
  expect_true(all(none$gene_table$resamples_qualified == 0))
  expect_true(all(is.na(none$gene_table$qualified_presence_fraction)))
  expect_equal(none$coverage_table$qualifying_resamples, 0)
  expect_true(is.na(none$gene_table$qualification_fraction[none$gene_table$gene_id == "g4"]))
})

test_that("lost universe metadata and inconsistent summaries are not guessed", {
  fixture <- resampling_fixture()
  attr(fixture$genes_table, "gene_order") <- NULL
  expect_error(summarize_gene_resampling(fixture), "complete gene_universe")
  explicit <- summarize_gene_resampling(fixture, gene_universe = c("g4", "g3", "g2", "g1"))
  expect_equal(explicit$coverage_table$sampling_coverage, c(.75, .75))
  fixture <- resampling_fixture()
  fixture$eval_table$summary_test_accuracy[2] <- .7
  expect_error(summarize_gene_resampling(fixture), "consistent summary")
  expect_error(summarize_gene_resampling(resampling_fixture(), algorithm = "unknown"), "Unknown")
  expect_error(summarize_gene_resampling(resampling_fixture(), gene_universe = c("g1", "g2")), "gene universe")
})
