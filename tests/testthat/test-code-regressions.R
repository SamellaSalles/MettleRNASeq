counts_fixture <- function(replicates = 4L, genes = 8L) {
  set.seed(42)
  x <- matrix(stats::rpois(2L * replicates * genes, 30), 2L * replicates, genes)
  x[(replicates + 1L):(2L * replicates), ] <- x[(replicates + 1L):(2L * replicates), ] + 20
  dimnames(x) <- list(paste0("sample", seq_len(nrow(x))), paste0("gene", seq_len(genes)))
  as.data.frame(x)
}
prepare_fixture <- function(n = 4L) {
  x <- counts_fixture(n)
  prepare_input_matrix(x, sample_vector = rownames(x), condition_vector = rep(c("A", "B"), each = n))
}

test_that("orientation uses sample labels and retains replicate identities", {
  x <- counts_fixture(4, 2)
  result <- prepare_input_matrix(t(x), sample_vector = rownames(x), condition_vector = rep(c("A", "B"), each = 4))
  expect_equal(dim(result), c(4L, 6L))
  expect_equal(result$Gene, rep(colnames(x), 2))
  expect_equal(attr(result, "sample_map")$Sample, rownames(x))
  expect_error(prepare_input_matrix(x, sample_vector = rownames(x)[-1], condition_vector = rep("A", 7)), "exactly one axis")
})

test_that("normalization keeps dimensions and names through filtering", {
  x <- data.frame(g1 = c(0, 10, 20, 30), g2 = c(0, 30, 20, 10), row.names = paste0("s", 1:4))
  result <- normalize_count_matrix(x, min_count = 1, group_size = 2, method = "scale")
  expect_identical(rownames(result), c("s2", "s3", "s4"))
  expect_identical(attr(result, "removed_samples"), "s1")
  one <- normalize_count_matrix(data.frame(g1 = 10:13, constant = rep(20, 4)), method = "scale")
  expect_equal(dim(one), c(4L, 1L))
  expect_equal(one$g1, as.numeric(scale(10:13, center = FALSE)))
  expect_error(normalize_count_matrix(x, group_size = 10), "exceeds")
})

test_that("TMM and DESeq2 preserve orientation on real engines", {
  x <- counts_fixture()
  for (method in c("TMM", "DESeq2")) {
    result <- normalize_count_matrix(x, min_count = 1, group_size = 2, method = method)
    expect_identical(rownames(result), rownames(x))
    expect_identical(colnames(result), colnames(x))
    expect_true(all(is.finite(as.matrix(result))))
  }
})

test_that("PCA selection matches identities and validates counts", {
  x <- counts_fixture()
  p <- prcomp(x)
  first <- feature_selection_pca(p, x, num_pcs = 2, num_genes_per_pc = 2)
  shuffled <- feature_selection_pca(p, x[,rev(seq_len(ncol(x)))], num_pcs = 2, num_genes_per_pc = 2)
  expect_equal(first, shuffled)
  expect_error(feature_selection_pca(p, x, num_pcs = 0), "integer")
})

test_that("classification works with 3, 4 and 5 replicates and one resample", {
  for (n in 3:5) {
    result <- classification_approach_advanced(prepare_fixture(n), "rf", list(rf = data.frame(mtry = 1)),
      num_resamples = 1, num_genes = 3, train_percentage = (n - 1)/n)
    expect_equal(dim(result$genes_table), c(3L, 1L))
    expect_equal(nrow(result$eval_table), n)
    expect_true(all(result$status_table$status == "ok"))
    expect_equal(nrow(result$test_eval_table), n)
    expect_equal(result$eval_table$summary_test_accuracy,
      rep(median(result$eval_table$test_accuracy), n))
  }
})

test_that("multiple held-out replicates are evaluated independently", {
  result <- classification_approach_advanced(prepare_fixture(), "rf", list(rf = data.frame(mtry = 1)),
    num_resamples = 1, num_genes = 3, train_percentage = .5, summary_stat = "mean")
  expect_equal(nrow(result$eval_table), 6L)
  expect_equal(nrow(result$test_eval_table), 12L)
  for (s in 1:6) expect_equal(result$eval_table$test_accuracy[s],
    mean(result$test_eval_table$test_accuracy[result$test_eval_table$data_split == s]))
  expect_equal(unique(result$eval_table$summary_test_accuracy), mean(result$eval_table$test_accuracy))
})

test_that("gene matching is robust to shuffled condition blocks", {
  input <- prepare_fixture()
  first <- classification_approach_advanced(input, "rf", list(rf = data.frame(mtry = 1)), 1, 3)
  # Keep the first condition gene order while scramble the other block.
  shuffled <- input[c(1:8,16:9), ]
  second <- classification_approach_advanced(shuffled, "rf", list(rf = data.frame(mtry = 1)), 1, 3)
  expect_equal(first$eval_table, second$eval_table)
  expect_equal(first$gene_ids_table, second$gene_ids_table)
})

test_that("bestTune matching supports non-mtry and multiple-parameter algorithms", {
  result <- classification_approach_advanced(prepare_fixture(), c("svmLinear2", "ranger"),
    list(svmLinear2 = data.frame(cost = c(.1,1)),
      ranger = expand.grid(mtry = 1, splitrule = "gini", min.node.size = c(1,2))), 1, 3)
  expect_equal(nrow(result$eval_table), 8L)
  expect_true(all(is.finite(result$eval_table$train_accuracy)))
})

test_that("failed model fits cannot silently produce successful summaries", {
  expect_error(classification_approach_advanced(prepare_fixture(), "not_a_model",
    list(not_a_model = data.frame(x = 1)), 1, 3), "resample 1, partition 1")
  expect_warning(result <- classification_approach_advanced(prepare_fixture(), "not_a_model",
    list(not_a_model = data.frame(x = 1)), 1, 3, error_handling = "continue"), "failed")
  expect_equal(nrow(result$eval_table), 4L)
  expect_true(all(result$status_table$status == "error"))
  expect_true(all(is.na(result$eval_table$summary_test_accuracy)))
})

test_that("gene selection uses IDs when matrix columns are reordered", {
  x <- counts_fixture()
  result <- classification_approach_advanced(prepare_fixture(), "rf", list(rf = data.frame(mtry = 1)), 1, 3)
  normal <- genes_from_classification_results(result$eval_table, result$genes_table, x, resamples = 1, return_plot = FALSE)
  reverse <- genes_from_classification_results(result$eval_table, result$genes_table, x[,8:1], resamples = 1, return_plot = FALSE)
  expect_equal(normal$filtered_matrix, reverse$filtered_matrix)
  expect_equal(normal$selected_gene_ids, unique(result$gene_ids_table[[1]]))
})

test_that("summary plots use recorded means and separate algorithms", {
  eval <- data.frame(algorithm = rep(c("rf", "svm"), each = 2), resample = 1,
    test_accuracy = c(.2,.8,.1,.3), summary_test_accuracy = c(.5,.5,.2,.2))
  genes <- data.frame(Resample_1 = 1:2)
  grDevices::pdf(file = tempfile(fileext = ".pdf")); on.exit(grDevices::dev.off())
  result <- genes_from_classification_results(eval, genes, counts_fixture(), resamples = 1)
  summary <- result$plot$layers[[1]]$data
  expect_equal(summary$summary_test_accuracy, c(.5,.2))
  expect_equal(nrow(summary), 2L)
})

test_that("complete gene sets and expression states have separate comparison modes", {
  first <- data.frame(LHS = c("{G1=low,G2=high}", "{G3=(0,10]}"), RHS = c("{G4=high}", "{G4=low}"), support = c(.5,.6))
  second <- data.frame(LHS = "{G2=low,G1=high}", RHS = "{G4=low}", support = .9)
  genes <- compare_gene_rules(first, second, plot = FALSE)
  expect_equal(genes$shared_rules$LHS, "{G1,G2}")
  expect_equal(genes$unique_rules_condition1$LHS, "{G3}")
  expect_setequal(genes$clusters$shared$Gene, c("G1","G2","G4"))
  expect_identical(genes$original_rules$condition1, first)
  items <- compare_gene_rules(first, second, comparison = "items", plot = FALSE)
  expect_equal(nrow(items$shared_rules), 0L)
  reordered <- data.frame(LHS = "{G2=high,G1=low}", RHS = "{G4=high}", support = .9)
  expect_equal(nrow(compare_gene_rules(first, reordered, comparison = "items", plot = FALSE)$shared_rules), 1L)
})

test_that("rule comparison handles empty inputs and interval commas", {
  expect_equal(nrow(compare_gene_rules(data.frame(),data.frame(),plot=FALSE)$shared_rules),0L)
  df <- data.frame(LHS = "{G1=(0,10],G2=(20,30]}", RHS = "{G3=high}")
  result <- compare_gene_rules(df, df, plot = FALSE, comparison = "items")
  expect_equal(result$shared_rules$LHS, df$LHS)
  expect_setequal(result$clusters$shared$Gene, c("G1","G2","G3"))
})

test_that("correlation differences align both axes and return unique pairs", {
  a <- matrix(c(1,.2,.2,1),2,dimnames=list(c("g1","g2"),c("g1","g2")))
  result <- icc_difference(list(correlation_matrix_cond1=a,correlation_matrix_cond2=a[,2:1]))
  expect_equal(result$absolute_diff_matrix, a*0)
  expect_equal(nrow(result$difference_table),1L)
  expect_equal(nrow(icc_difference(list(correlation_matrix_cond1=a,correlation_matrix_cond2=a), unique_pairs=FALSE)$difference_table),4L)
  expect_error(acc_analysis(counts_fixture(),c(3,8)), "equal")
  expect_error(icc_analysis(counts_fixture(),c(0,8)), "two groups")
})

test_that("target-gene filters do not match partial gene names", {
  x <- data.frame(G1 = factor(c("low","low","high","high")),
    G10 = factor(c("low","low","high","high")), G2 = factor(c("low","low","high","high")))
  grDevices::pdf(file = tempfile(fileext = ".pdf")); on.exit(grDevices::dev.off())
  result <- mine_gene_associations(x,1:4,supp=.5,conf=.8,maxlen=2,target_genes="G1",search_type="lhs",interactive=FALSE)
  sides <- as(arules::lhs(result$rules),"list")
  expect_gt(length(sides),0L)
  expect_false(is.null(result$plot))
  expect_true(all(vapply(sides,function(items) any(sub("=.*$","",items)=="G1"),logical(1))))
})

test_that("one-gene mining and empty plots retain valid structures", {
  x <- data.frame(G1=factor(c("low","low","high","high")))
  expect_equal(length(mine_gene_associations(x,1:4,maxlen=2,interactive=FALSE)$rules),0L)
  grDevices::pdf(file = tempfile(fileext = ".pdf")); on.exit(grDevices::dev.off())
  expect_silent(compare_gene_rules(data.frame(),data.frame(),plot=TRUE))
})

test_that("PCA plots support factor metadata and numeric colors", {
  grDevices::pdf(file = tempfile(fileext = ".pdf")); on.exit(grDevices::dev.off())
  x <- counts_fixture()
  expect_no_error(pca_analysis(x, color_by = factor(rep(c("A", "B"), each = 4)),
    shape_by = rep(1:2, each = 4), plot_3d = TRUE))
  expect_no_error(pca_analysis(x, color_by = seq_len(nrow(x))))
  expect_no_error(pca_analysis(x, plot_3d = TRUE))
})

test_that("rule graphs can draw more than twelve disconnected clusters", {
  grDevices::pdf(file = tempfile(fileext = ".pdf")); on.exit(grDevices::dev.off())
  rules <- data.frame(LHS = paste0("{L", 1:15, "=low}"), RHS = paste0("{R", 1:15, "=high}"))
  expect_no_error(result <- compare_gene_rules(rules, rules, plot = TRUE))
  expect_equal(nrow(result$clusters$shared), 30L)
  expect_equal(length(unique(result$clusters$shared$Cluster)), 15L)
})

test_that("FP-growth requires an explicit optional backend and valid thresholds", {
  x <- data.frame(G1 = factor(c("low", "low", "high", "high")),
    G2 = factor(c("low", "low", "high", "high")))
  if (!requireNamespace("fim4r", quietly = TRUE)) {
    expect_error(mine_gene_associations(x, 1:4, method = "fpgrowth"), "optional fim4r")
  } else {
    grDevices::pdf(file = tempfile(fileext = ".pdf")); on.exit(grDevices::dev.off())
    result <- mine_gene_associations(x, 1:4, method = "fpgrowth", supp = .5,
      conf = .9, maxlen = 2, interactive = FALSE)
    expect_gt(length(result$rules), 0L)
    expect_true(all(arules::quality(result$rules)$support >= .5))
    expect_true(all(arules::quality(result$rules)$confidence >= .9))
  }
})
