  #' Clean and Normalize a Count Matrix Using TMM, DESeq2, or Scaling
  #'
  #' This function normalizes a count matrix using one of three methods: TMM (edgeR), DESeq2, or scaling.
  #' It also provides options to filter out low-variance genes, remove low-count genes, and remove samples with zero library size.
  #'
  #' @param count_matrix A data frame or matrix where rows are samples and columns are genes. Row names should be sample names.
  #' @param remove_low_variance Logical. If `TRUE`, genes with zero variance across samples are removed. Default is `TRUE`.
  #' @param min_count Numeric. The minimum count threshold for a gene to be retained. Genes must meet this threshold in at least `group_size` samples to be retained. Default is `10`.
  #' @param group_size Numeric. The minimum number of samples where a gene must meet the `min_count` threshold to be retained. Default is `4`.
  #' @param method Character. The normalization method to use. Options are `"TMM"` (edgeR), `"DESeq2"`, or `"scale"`. Default is `"TMM"`.
  #'
  #' @return A normalized count matrix (data frame) with samples as rows and genes as columns. Row names are preserved.
  #'
  #' @details
  #' The function performs the following steps:
  #' 1. Converts the input matrix to numeric and preserves row names.
  #' 2. Removes genes with zero variance if `remove_low_variance = TRUE`.
  #' 3. Filters out genes with counts below `min_count` in fewer than `group_size` samples.
  #' 4. Removes samples with zero library size.
  #' 5. Normalizes the count matrix using the specified method:
  #'    - **TMM**: Trimmed Mean of M-values normalization (edgeR).
  #'    - **DESeq2**: Size factor normalization (DESeq2).
  #'    - **scale**: Non-centered column-wise scaling, preserving the original scale(center = FALSE, scale = TRUE) convention.
  #'
  #' @examples
  #' # Example count matrix (samples as rows, genes as columns)
  #' count_matrix <- data.frame(
  #'   Gene1 = c(10, 20, 30, 40),
  #'   Gene2 = c(15, 25, 35, 45),
  #'   Gene3 = c(12, 22, 32, 42)
  #' )
  #'
  #' # Normalize using TMM
  #' norm_counts <- normalize_count_matrix(count_matrix, method = "TMM")
  #' 
  #' # Normalize using DESeq2
  #' norm_counts <- normalize_count_matrix(count_matrix, method = "DESeq2")
  #'
  #' # Normalize using scaling
  #' norm_counts <- normalize_count_matrix(count_matrix, method = "scale")
  #'
  #' @importFrom edgeR DGEList calcNormFactors
  #' @importFrom DESeq2 DESeqDataSetFromMatrix estimateSizeFactors counts
  #' @importFrom dplyr %>%
  #' @export
  normalize_count_matrix <- function(count_matrix, remove_low_variance = TRUE, min_count = 10,
                                   group_size = 4, method = "TMM") {
  method <- match.arg(method, c("TMM", "DESeq2", "scale"))
  group_size <- .mettle_integer(group_size, "group_size")
  if (length(min_count) != 1L || !is.numeric(min_count) || !is.finite(min_count) || min_count < 0) stop("min_count must be finite and non-negative.")
  if (!is.logical(remove_low_variance) || length(remove_low_variance) != 1L || is.na(remove_low_variance)) stop("remove_low_variance must be TRUE or FALSE.")
  x <- .mettle_numeric_matrix(count_matrix, "count_matrix", nonnegative = TRUE)
  if (group_size > nrow(x)) stop("group_size exceeds the number of samples.")
  if (is.null(rownames(x))) rownames(x) <- paste0("Sample_", seq_len(nrow(x)))
  if (is.null(colnames(x))) colnames(x) <- paste0("Gene_", seq_len(ncol(x)))
  .mettle_names(rownames(x), "Sample names"); .mettle_names(colnames(x), "Gene names")
  original_samples <- rownames(x); original_genes <- colnames(x)
  keep <- colSums(x >= min_count) >= group_size
  if (remove_low_variance) {
    if (nrow(x) < 2L) stop("Variance filtering requires at least two samples.")
    keep <- keep & apply(x, 2, stats::var) > 0
  }
  x <- x[, keep, drop = FALSE]
  x <- x[rowSums(x) > 0, , drop = FALSE]
  if (!nrow(x) || !ncol(x)) stop("No genes or samples left after filtering. Adjust min_count or group_size.")
  if (method == "TMM") {
    dge <- edgeR::calcNormFactors(edgeR::DGEList(counts = t(x)), method = "TMM")
    effective_sizes <- dge$samples$lib.size * dge$samples$norm.factors
    size_factors <- effective_sizes / exp(mean(log(effective_sizes)))
    result <- t(sweep(t(x), 2, size_factors, "/"))
  } else if (method == "DESeq2") {
    if (any(x != round(x))) stop("DESeq2 normalization requires integer counts.")
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = t(x),
      colData = data.frame(row.names = rownames(x)), design = ~ 1)
    result <- t(DESeq2::counts(DESeq2::estimateSizeFactors(dds), normalized = TRUE))
  } else {
    # Preserve the original non-centered scaling convention.
    result <- base::scale(x, center = FALSE, scale = TRUE)
  }
  out <- as.data.frame(result)
  attr(out, "removed_samples") <- setdiff(original_samples, rownames(x))
  attr(out, "removed_genes") <- setdiff(original_genes, colnames(x))
  out
}
