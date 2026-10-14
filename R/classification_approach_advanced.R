#' Advanced Classification with Gene Resampling and Replicate Splitting
#'
#' This function performs classification using selected algorithms across multiple resamples of genes,
#' with all possible train/test replicate partitions. It outputs evaluation metrics, selected genes,
#' partitioning details, and trained models.
#'
#' @param input_data A data frame output from \code{prepare_input_matrix()}, with genes as rows,
#' replicates (e.g., \code{Rep_1}, \code{Rep_2}, ...) as columns, and a \code{Condition} column.
#' @param algorithms A character vector of classification algorithms (e.g., \code{"rf"}, \code{"ranger"}, \code{"svmLinear2"}).
#' @param hyper_params A named list of tuning grids (data.frames) for each algorithm. Each list element must contain
#' the hyperparameters specific to the corresponding algorithm. See details below.
#' @param num_resamples Number of bootstrap resamples (i.e., how many different sets of genes to sample).
#' @param num_genes Number of genes to randomly sample from the first condition in each resample.
#' @param train_percentage Percentage of replicates to use for training (e.g., 0.75).
#' @param metrics A character vector of evaluation metrics (currently only \code{"Accuracy"} is supported).
#' @param summary_stat Statistic to summarize test performance per resample. Either \code{"median"} or \code{"mean"}.
#' @param replace_genes Logical; whether to sample genes with replacement.
#' @param parallel Logical; whether to allow parallel execution during model training (if supported by the algorithm).
#' @param seed Integer; random seed for reproducibility.
#'
#' @param include_gene Include gene identity as a predictor (TRUE preserves the original approach).
#' @param error_handling Stop at a failed fit (default), or continue with explicit failure records and NA summaries.
#' @details
#' Genes are sampled only from the first condition, then matched across all conditions. For each resample, 
#' the function uses all possible combinations of train/test replicate splits based on the given training percentage.
#'
#' The \code{hyper_params} argument must be a named list where each name corresponds to an algorithm (e.g., \code{"rf"} or \code{"ranger"}), and each value is a data frame with hyperparameter settings. Example:
#'
#' \preformatted{
#' hyper_params <- list(
#'   rf = data.frame(mtry = c(1, 5, 10, sqrt(length(cols)))), 
#'   ranger = data.frame(mtry = c(1, 5, 10, sqrt(length(cols))), 
#'                       splitrule = "extratrees", 
#'                       min.node.size = c(2, 3))
#' )
#' }
#'
#' Each held-out replicate is evaluated separately by repeating its expression values across the training replicate columns. Partition accuracy is the mean of these replicate accuracies. No internal gene holdout is introduced. Gene identity remains a predictor by default; a constant identifier is omitted. Additional outputs are gene_ids_table, status_table, test_eval_table, sample_map and summary_stat.
#'
#' Output metrics include training and test accuracy per split and per resample, with summaries per resample.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{eval_table}{A data frame with algorithm name, resample number, data split, train/test replicates, and metrics.}
#'   \item{partition_table}{A data frame listing all unique train/test replicate combinations used.}
#'   \item{genes_table}{A data frame where each column is a resample and each row is the position (index) of a gene selected from the first condition.}
#'   \item{models}{A nested list of trained models per algorithm, resample, and data split.}
#' }
#'
#' @examples
#' \dontrun{
#' # Example sample and condition vectors
#' sample <- c("F1B", "F2B", "F3B", "F4B", "S1B", "S2B", "S3B", "S4B")
#' condition <- c("FLASH", "FLASH", "FLASH", "FLASH", "STAN", "STAN", "STAN", "STAN")
#' data <- matrix(runif(8 * 100), nrow = 100, ncol = 8)
#' colnames(data) <- sample
#' result <- prepare_input_matrix(count_matrix = data, sample_vector = sample, condition_vector = condition)
#' result
#' 
#' hyper_params <- list(
#'   rf = data.frame(mtry = c(1, 5, 10))
#' )
#' results <- classification_approach_advanced(
#'   input_data = result,
#'   algorithms = c("rf"),
#'   hyper_params = hyper_params,
#'   num_resamples = 5,
#'   num_genes = 50
#' )
#' }
#'
#' @export
classification_approach_advanced <- function(input_data, algorithms, hyper_params, num_resamples, num_genes,
                                             train_percentage = 0.75, metrics = "Accuracy", summary_stat = "median",
                                             replace_genes = TRUE, parallel = FALSE, seed = 123,
                                             include_gene = TRUE, error_handling = c("stop", "continue")) {
  error_handling <- match.arg(error_handling)
  summary_stat <- match.arg(summary_stat, c("median", "mean"))
  num_resamples <- .mettle_integer(num_resamples, "num_resamples")
  num_genes <- .mettle_integer(num_genes, "num_genes")
  seed <- .mettle_integer(seed, "seed", 0L)
  if (!identical(metrics, "Accuracy")) stop("Only metrics = 'Accuracy' is supported.")
  for (flag in list(replace_genes, parallel, include_gene)) {
    if (!is.logical(flag) || length(flag) != 1L || is.na(flag)) stop("Logical flags must be TRUE or FALSE.")
  }
  if (!is.data.frame(input_data) || !all(c("Gene", "Condition") %in% names(input_data))) {
    stop("input_data must contain Gene and Condition columns.")
  }
  .mettle_names(names(input_data), "Input column names")
  if (anyNA(input_data$Gene) || anyNA(input_data$Condition) ||
      any(!nzchar(as.character(input_data$Gene))) || any(!nzchar(as.character(input_data$Condition)))) {
    stop("Gene and Condition must contain non-empty labels without NA.")
  }
  replicate_cols <- setdiff(names(input_data), c("Gene", "Condition"))
  if (length(replicate_cols) < 2L) stop("At least two replicate columns are required.")
  .mettle_numeric_matrix(input_data[, replicate_cols, drop = FALSE], "Replicate values")
  conditions <- unique(as.character(input_data$Condition))
  if (length(conditions) < 2L) stop("At least two conditions are required.")
  blocks <- lapply(conditions, function(c) which(as.character(input_data$Condition) == c))
  gene_ids <- as.character(input_data$Gene[blocks[[1]]])
  .mettle_names(gene_ids, "Genes within each condition")
  # Match by identity rather than assuming contiguous, identically ordered blocks.
  blocks <- lapply(blocks, function(idx) {
    ids <- as.character(input_data$Gene[idx])
    .mettle_names(ids, "Genes within each condition")
    if (!setequal(ids, gene_ids)) stop("All conditions must contain exactly the same genes.")
    idx[match(gene_ids, ids)]
  })
  if (!replace_genes && num_genes > length(gene_ids)) stop("num_genes exceeds available genes without replacement.")
  if (!is.numeric(train_percentage) || length(train_percentage) != 1L ||
      !is.finite(train_percentage) || train_percentage <= 0 || train_percentage >= 1) {
    stop("train_percentage must be strictly between 0 and 1.")
  }
  n_train <- round(train_percentage * length(replicate_cols))
  if (n_train < 1L || n_train >= length(replicate_cols)) stop("The requested split must leave both training and test replicates.")
  if (!is.character(algorithms) || !length(algorithms) || anyNA(algorithms) ||
      any(!nzchar(algorithms)) || anyDuplicated(algorithms)) stop("algorithms must contain unique method names.")
  if (!is.list(hyper_params) || !all(algorithms %in% names(hyper_params))) stop("Supply a tuning grid for every algorithm.")
  for (m in algorithms) {
    if (!is.data.frame(hyper_params[[m]]) || !nrow(hyper_params[[m]])) stop("Each tuning grid must be a non-empty data frame.")
  }
  if (!requireNamespace("caret", quietly = TRUE)) stop("Install caret before classification.")
  splits <- utils::combn(seq_along(replicate_cols), n_train, simplify = FALSE)
  partition_table <- do.call(rbind, lapply(seq_along(splits), function(s) data.frame(
    Partition_ID = paste("Partition", s),
    Train_Splits = paste(replicate_cols[splits[[s]]], collapse = ","),
    Test_Splits = paste(replicate_cols[setdiff(seq_along(replicate_cols), splits[[s]])], collapse = ","))))
  set.seed(seed)
  selections <- lapply(seq_len(num_resamples), function(i) sample(seq_along(gene_ids), num_genes, replace = replace_genes))
  genes_table <- as.data.frame(do.call(cbind, selections))
  names(genes_table) <- paste0("Resample_", seq_len(num_resamples))
  gene_ids_table <- as.data.frame(lapply(selections, function(idx) gene_ids[idx]), stringsAsFactors = FALSE)
  names(gene_ids_table) <- names(genes_table)
  attr(genes_table, "gene_ids") <- gene_ids_table
  attr(genes_table, "gene_order") <- gene_ids
  eval_rows <- list(); test_rows <- list(); status_rows <- list(); models <- list()
  for (m in algorithms) {
    models[[m]] <- list()
    for (i in seq_len(num_resamples)) {
      idx <- unlist(lapply(blocks, function(b) b[selections[[i]]]), use.names = FALSE)
      subset <- input_data[idx, , drop = FALSE]
      reference <- factor(as.character(subset$Condition), levels = conditions)
      for (s in seq_along(splits)) {
        train_cols <- replicate_cols[splits[[s]]]
        test_cols <- setdiff(replicate_cols, train_cols)
        key <- paste0("resample_", i, "_split_", s)
        trained <- tryCatch({
          training <- subset[, train_cols, drop = FALSE]
          # A constant identifier has no information and cannot be dummy-encoded.
          use_gene <- include_gene && length(unique(subset$Gene)) > 1L
          if (use_gene) training <- cbind(Gene = as.character(subset$Gene), training)
          training$Condition <- reference
          control <- caret::trainControl(method = "boot", number = 1,
            index = list(seq_len(nrow(training))), indexOut = list(seq_len(nrow(training))),
            allowParallel = parallel)
          set.seed(seed)
          model <- caret::train(Condition ~ ., data = training, method = m, metric = "Accuracy",
            preProcess = "center", trControl = control, tuneGrid = hyper_params[[m]])
          # Match every tuning parameter, not only mtry.
          best <- rep(TRUE, nrow(model$results))
          for (param in names(model$bestTune)) {
            value <- model$bestTune[[param]][1]
            best <- best & if (is.na(value)) is.na(model$results[[param]]) else
              !is.na(model$results[[param]]) & model$results[[param]] == value
          }
          accuracy <- model$results$Accuracy[best]
          if (length(accuracy) != 1L || !is.finite(accuracy)) stop("Could not identify one finite training accuracy for bestTune.")
          per_test <- vapply(test_cols, function(test_col) {
            # Evaluate each held-out replicate separately, repeating its values as
            # in the default leave-one-replicate-out procedure.
            testing <- as.data.frame(matrix(rep(subset[[test_col]], length(train_cols)), ncol = length(train_cols)))
            names(testing) <- train_cols
            if (use_gene) testing <- cbind(Gene = as.character(subset$Gene), testing)
            predicted <- stats::predict(model, newdata = testing)
            if (length(predicted) != length(reference) || anyNA(predicted)) stop("Invalid predictions.")
            mean(as.character(predicted) == as.character(reference))
          }, numeric(1))
          list(model = model, train_accuracy = accuracy, per_test = per_test)
        }, error = function(e) e)
        ok <- !inherits(trained, "error")
        error <- if (ok) NA_character_ else conditionMessage(trained)
        status_rows[[length(status_rows) + 1L]] <- data.frame(algorithm = m, resample = i,
          data_split = s, status = if (ok) "ok" else "error", error = error)
        if (!ok && error_handling == "stop") stop(sprintf("%s, resample %d, partition %d: %s", m, i, s, error), call. = FALSE)
        eval_rows[[length(eval_rows) + 1L]] <- data.frame(algorithm = m, resample = i,
          data_split = s, train_partition = paste(train_cols, collapse = ","),
          test_partition = paste(test_cols, collapse = ","),
          train_accuracy = if (ok) trained$train_accuracy else NA_real_,
          test_accuracy = if (ok) mean(trained$per_test) else NA_real_,
          status = if (ok) "ok" else "error", error = error)
        if (ok) {
          models[[m]][[key]] <- trained$model
          test_rows[[length(test_rows) + 1L]] <- data.frame(algorithm = m, resample = i, data_split = s,
            test_replicate = test_cols, test_accuracy = unname(trained$per_test))
        }
      }
    }
  }
  eval_table <- do.call(rbind, eval_rows)
  eval_table$summary_test_accuracy <- NA_real_
  for (m in algorithms) for (i in seq_len(num_resamples)) {
    rows <- which(eval_table$algorithm == m & eval_table$resample == i)
    # Never present a failed/incomplete set of partitions as a complete score.
    if (all(eval_table$status[rows] == "ok")) {
      eval_table$summary_test_accuracy[rows] <- if (summary_stat == "median")
        stats::median(eval_table$test_accuracy[rows]) else mean(eval_table$test_accuracy[rows])
    }
  }
  status_table <- do.call(rbind, status_rows)
  if (any(status_table$status == "error")) warning("Some model fits failed; see status_table. Incomplete resamples have NA summaries.")
  list(eval_table = eval_table, partition_table = partition_table, genes_table = genes_table,
    models = models, gene_ids_table = gene_ids_table, status_table = status_table,
    test_eval_table = if (length(test_rows)) do.call(rbind, test_rows) else data.frame(),
    sample_map = attr(input_data, "sample_map"), summary_stat = summary_stat)
}
