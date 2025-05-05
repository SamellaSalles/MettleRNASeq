  #' Clean and Normalize a Count Matrix Using TMM, DESeq2, or Scaling
  #'
  #' This function normalizes a count matrix using one of three methods: TMM (edgeR), DESeq2, or scaling.
  #' It also provides options to filter out low-variance genes, remove low-count genes, and remove samples with zero library size.
  #'
  #' @param count_matrix A data frame or matrix where rows are samples and columns are genes. Row names should be sample names.
  #' @param remove_low_variance Logical. If `TRUE`, genes with zero variance across samples are removed. Default is `TRUE`.
  #' @param min_count Numeric. The minimum count threshold for a gene to be retained. Genes with counts below this threshold in at least `group_size` samples are removed. Default is `10`.
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
  #'    - **scale**: Column-wise scaling (mean-centered and scaled to unit variance).
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
  normalize_count_matrix <- function(count_matrix, remove_low_variance = TRUE, min_count = 10, group_size = 4, method = "TMM") {
    require(edgeR)
    require(DESeq2)
    require(dplyr)
    
    # Ensure samples are in row names and genes are in columns
    if (!is.data.frame(count_matrix)) {
      stop("Error: count_matrix must be a data frame with samples as row names and genes as columns.")
    }
    
    # Preserve original row names (samples)
    sample_names <- rownames(count_matrix)
    
    # Convert to numeric if necessary
    count_matrix <- as.data.frame(lapply(count_matrix, as.numeric))
    rownames(count_matrix) <- sample_names  # Restore row names
    
    # Remove genes with low variance if selected
    if (remove_low_variance) {
      keep_genes <- apply(count_matrix, 2, function(x) var(x, na.rm = TRUE) > 0)
      count_matrix <- count_matrix[, keep_genes]
    }
    
    # Allow user to specify minimum count and group size
    keep <- colSums(count_matrix >= min_count) >= group_size
    count_matrix <- count_matrix[, keep]
    
    # Remove samples (rows) with zero library size
    non_zero_samples <- rowSums(count_matrix) > 0
    count_matrix <- count_matrix[non_zero_samples, ]
    
    # Check if matrix is empty after filtering
    if (nrow(count_matrix) == 0 || ncol(count_matrix) == 0) {
      stop("Error: No genes or samples left after filtering. Consider adjusting min_count or group_size.")
    }
    
    # Normalization methods
    if (method == "TMM") {
      dge <- DGEList(counts = t(count_matrix))  # Transpose to match edgeR input format
      dge <- calcNormFactors(dge, method = "TMM")
      
      norm_factors <- dge$samples$norm.factors
      lib_sizes <- colSums(t(count_matrix))
      size_factors <- lib_sizes * norm_factors / exp(mean(log(lib_sizes * norm_factors)))
      
      norm_counts <- sweep(t(count_matrix), 2, size_factors, "/")  # Apply normalization
      norm_counts <- t(norm_counts)  # Transpose back to original format
    } else if (method == "DESeq2") {
      dds <- DESeqDataSetFromMatrix(countData = t(count_matrix), colData = data.frame(row.names = rownames(count_matrix)), design = ~ 1)
      dds <- estimateSizeFactors(dds)
      norm_counts <- counts(dds, normalized = TRUE)
      norm_counts <- t(norm_counts)  # Transpose back to original format
    } else if (method == "scale") {
      norm_counts <- scale(count_matrix, center = FALSE, scale = TRUE)
    } else {
      stop("Invalid normalization method. Choose 'TMM', 'DESeq2', or 'scale'.")
    }
    
    # Restore original row names
    rownames(norm_counts) <- sample_names
    
    # Convert to data frame
    norm_counts <- as.data.frame(norm_counts)
    
    return(norm_counts)
  }
  