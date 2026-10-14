#' Prepare Input Matrix for Classification
#'
#' This function prepares an input matrix for classification by organizing count data
#' into a wide format with consistent replicate naming.
#'
#' @param count_matrix A numeric matrix or data frame with genes on one axis and samples on the other. Sample labels must match exactly one axis; orientation is determined from those labels.
#' @param design_table A data frame containing sample and condition information. Default is NULL.
#' @param sample_col The name of the column in `design_table` containing sample names. Default is NULL.
#' @param condition_col The name of the column in `design_table` containing condition names. Default is NULL.
#' @param sample_vector A vector of sample names. Default is NULL.
#' @param condition_vector A vector of condition names corresponding to `sample_vector`. Default is NULL.
#'
#' @return A data frame with one row per gene and condition, a Gene column,
#' replicate columns named Rep_1, Rep_2, etc., and a final Condition column.
#' Its sample_map attribute links each condition/replicate slot to its original sample ID.
#'
#' @details
#' This function:
#' \itemize{
#'   \item Ensures all conditions have the same number of replicates
#'   \item Converts the input matrix to a properly formatted data frame
#'   \item Uses underscore-separated replicate names ("Rep_1" instead of "Replicate 1")
#'   \item Maintains original gene order while arranging by condition
#' }
#'
#' @examples
#' # Example with vectors
#' sample <- c("F1", "F2", "F3", "F4", "S1", "S2", "S3", "S4")
#' condition <- c("F", "F", "F", "F", "S", "S", "S", "S")
#' result <- matrix(runif(8 * 100), nrow = 100, ncol = 8)
#' colnames(result) <- sample
#' prepare_input_matrix(count_matrix = result, sample_vector = sample, condition_vector = condition)
#'
#' @importFrom dplyr select arrange group_by ungroup mutate relocate
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom stats setNames
#' @export
prepare_input_matrix <- function(count_matrix, design_table = NULL, sample_col = NULL,
                                 condition_col = NULL, sample_vector = NULL, condition_vector = NULL) {
  if (!is.null(design_table)) {
    if (is.null(sample_col) || is.null(condition_col) ||
        !all(c(sample_col, condition_col) %in% names(design_table))) {
      stop("Supply existing sample_col and condition_col columns in design_table.")
    }
    samples <- as.character(design_table[[sample_col]])
    conditions <- as.character(design_table[[condition_col]])
  } else {
    if (is.null(sample_vector) || is.null(condition_vector)) stop("Supply design_table or sample_vector and condition_vector.")
    samples <- as.character(sample_vector); conditions <- as.character(condition_vector)
  }
  .mettle_names(samples, "Sample labels")
  if (length(samples) != length(conditions) || anyNA(conditions) || any(!nzchar(conditions))) stop("Each sample must have one non-empty condition.")
  x <- .mettle_numeric_matrix(count_matrix, "count_matrix")
  in_cols <- !is.null(colnames(x)) && setequal(samples, colnames(x))
  in_rows <- !is.null(rownames(x)) && setequal(samples, rownames(x))
  if (in_cols && in_rows) stop("Ambiguous orientation: sample labels match both axes.")
  if (!in_cols && !in_rows) stop("Sample labels must match exactly one axis of count_matrix, without extra or missing samples.")
  if (in_rows) x <- t(x)
  .mettle_names(colnames(x), "Sample names")
  if (is.null(rownames(x))) rownames(x) <- paste0("Gene_", seq_len(nrow(x)))
  .mettle_names(rownames(x), "Gene names")
  # Preserve original matrix sample order, including the original replicate slots.
  conditions <- conditions[match(colnames(x), samples)]
  levels <- if (!is.null(design_table)) unique(as.character(design_table[[condition_col]])) else unique(as.character(condition_vector))
  counts <- table(factor(conditions, levels = levels))
  if (length(levels) < 2L || length(unique(as.integer(counts))) != 1L) stop("At least two conditions with equal replicate counts are required.")
  blocks <- lapply(levels, function(c) {
    idx <- which(conditions == c)
    values <- as.data.frame(x[, idx, drop = FALSE])
    names(values) <- paste0("Rep_", seq_along(idx))
    data.frame(Gene = rownames(x), values, Condition = c, check.names = FALSE)
  })
  out <- do.call(rbind, blocks); rownames(out) <- NULL
  attr(out, "sample_map") <- do.call(rbind, lapply(levels, function(c) {
    idx <- which(conditions == c)
    data.frame(Condition = c, Replicate = paste0("Rep_", seq_along(idx)), Sample = colnames(x)[idx])
  }))
  out
}
