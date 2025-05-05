#' Perform Feature Selection Using PCA Results
#'
#' This function performs feature selection based on the results of Principal Component Analysis (PCA).
#' It allows the user to select the number of principal components (PCs) and the number of top or worst genes to extract from each PC.
#' The function also allows the user to specify whether to select genes from the first PCs or the last PCs.
#' The function returns a new count matrix containing only the selected genes.
#'
#' @param pca_result The result of the `prcomp` function (PCA result).
#' @param count_matrix A data frame or matrix where rows are samples and columns are genes, or vice versa. Row names should be sample or gene names.
#' @param num_pcs Numeric. The number of principal components (PCs) to consider for feature selection. Default is 3.
#' @param num_genes_per_pc Numeric. The number of top or worst genes to extract from each selected PC. Default is 100.
#' @param select_worst Logical. If `TRUE`, select the worst genes (smallest absolute loadings). If `FALSE`, select the top genes (largest absolute loadings). Default is `FALSE`.
#' @param use_last_pcs Logical. If `TRUE`, select genes from the last `num_pcs` PCs. If `FALSE`, select genes from the first `num_pcs` PCs. Default is `FALSE`.
#'
#' @return A list containing:
#' \itemize{
#'   \item `selected_genes`: A vector of gene names selected based on the PCA results.
#'   \item `selected_count_matrix`: A new count matrix containing only the selected genes.
#' }
#'
#' @examples
#' # Example count matrix (samples as rows, genes as columns)
#' count_matrix <- data.frame(
#'   Gene1 = c(10, 20, 30, 40),
#'   Gene2 = c(15, 25, 35, 45),
#'   Gene3 = c(12, 22, 32, 42)
#' )
#'
#' # Perform PCA
#' pca_result <- prcomp(count_matrix, center = TRUE, scale = TRUE)
#'
#' # Perform feature selection (top genes from first PCs)
#' feature_selection_results <- feature_selection_pca(pca_result, count_matrix, num_pcs = 2, num_genes_per_pc = 1)
#'
#' # Access selected genes and count matrix
#' feature_selection_results$selected_genes
#' feature_selection_results$selected_count_matrix
#'
#' # Perform feature selection (worst genes from last PCs)
#' feature_selection_results_worst <- feature_selection_pca(pca_result, count_matrix, num_pcs = 2, num_genes_per_pc = 1, select_worst = TRUE, use_last_pcs = TRUE)
#'
#' # Access selected genes and count matrix
#' feature_selection_results_worst$selected_genes
#' feature_selection_results_worst$selected_count_matrix
#'
#' @importFrom stats prcomp
#' @export
feature_selection_pca <- function(pca_result, count_matrix, num_pcs = 3, num_genes_per_pc = 100, select_worst = FALSE, use_last_pcs = FALSE) {
  # Convert count_matrix to a data frame if it isn't already
  if (!is.data.frame(count_matrix)) {
    count_matrix <- as.data.frame(count_matrix)
  }
  
  # Ensure row names are preserved
  if (is.null(rownames(count_matrix))) {
    stop("Error: count_matrix must have row names.")
  }
  
  # Ensure pca_result is a valid PCA result
  if (!inherits(pca_result, "prcomp")) {
    stop("Error: pca_result must be the result of the prcomp function.")
  }
  
  # Ensure num_pcs and num_genes_per_pc are valid
  if (num_pcs > ncol(pca_result$rotation)) {
    stop("Error: num_pcs cannot be greater than the number of PCs in pca_result.")
  }
  if (num_genes_per_pc > nrow(pca_result$rotation)) {
    stop("Error: num_genes_per_pc cannot be greater than the number of genes in pca_result.")
  }
  
  # Check if count_matrix needs to be transposed
  # If the number of rows in count_matrix matches the number of genes in pca_result, assume genes are rows
  if (nrow(count_matrix) == nrow(pca_result$rotation)) {
    count_matrix <- t(count_matrix)  # Transpose to make samples rows and genes columns
  } else if (ncol(count_matrix) == nrow(pca_result$rotation)) {
    # Already in the correct orientation (samples as rows, genes as columns)
  } else {
    stop("Error: count_matrix dimensions do not match the PCA result. Ensure samples are rows and genes are columns.")
  }
  
  # Extract the rotation matrix (loadings) from the PCA result
  rotation_matrix <- pca_result$rotation
  
  # Determine which PCs to use
  if (use_last_pcs) {
    # Use the last `num_pcs` PCs
    pc_indices <- (ncol(rotation_matrix) - num_pcs + 1):ncol(rotation_matrix)
  } else {
    # Use the first `num_pcs` PCs
    pc_indices <- 1:num_pcs
  }
  
  # Initialize a vector to store selected genes
  selected_genes <- c()
  
  # Loop through the specified PCs
  for (i in pc_indices) {
    # Get the absolute loadings for the current PC
    pc_loadings <- abs(rotation_matrix[, i])
    
    # Sort the genes by their absolute loadings
    if (select_worst) {
      # Sort in ascending order for worst genes (smallest absolute loadings)
      sorted_genes <- names(sort(pc_loadings, decreasing = FALSE))
    } else {
      # Sort in descending order for top genes (largest absolute loadings)
      sorted_genes <- names(sort(pc_loadings, decreasing = TRUE))
    }
    
    # Select the top or worst genes for the current PC
    top_genes <- sorted_genes[1:num_genes_per_pc]
    
    # Add the selected genes to the list
    selected_genes <- c(selected_genes, top_genes)
  }
  
  # Remove duplicate genes (if any)
  selected_genes <- unique(selected_genes)
  
  # Create a new count matrix with only the selected genes
  selected_count_matrix <- count_matrix[, selected_genes, drop = FALSE]
  
  # Return the results
  return(list(
    selected_genes = selected_genes,
    selected_count_matrix = selected_count_matrix
  ))
}
