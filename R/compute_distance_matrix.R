#' Compute a Distance Matrix and Generate a Heatmap
#'
#' This function computes a distance matrix from a count matrix and optionally generates a heatmap.
#' The user can choose between distance-based methods (e.g., Euclidean, Manhattan) or correlation-based methods (e.g., Spearman, Pearson).
#' The heatmap can be customized with clustering, color palettes, and other parameters.
#'
#' @param count_matrix A data frame or matrix where rows are samples and columns are genes. Row names should be sample names.
#' @param use_correlation Logical. If `TRUE`, the distance matrix is computed using correlation. If `FALSE`, a distance-based method is used. Default is `FALSE`.
#' @param method Character. The distance method to use when `use_correlation = FALSE`. Options are `"euclidean"`, `"manhattan"`, `"canberra"`, or `"minkowski"`. Default is `"euclidean"`.
#' @param correlation_method Character. The correlation method to use when `use_correlation = TRUE`. Options are `"spearman"`, `"pearson"`, or `"kendall"`. Default is `"spearman"`.
#' @param plot_heatmap Logical. If `TRUE`, a heatmap of the distance matrix is generated. Default is `TRUE`.
#' @param cluster_rows Logical. If `TRUE`, rows of the heatmap are clustered. Default is `TRUE`.
#' @param cluster_cols Logical. If `TRUE`, columns of the heatmap are clustered. Default is `TRUE`.
#' @param clustering_method Character. The clustering method to use for rows and columns. Options include `"complete"`, `"average"`, `"single"`, `"ward.D"`, `"ward.D2"`, `"mcquitty"`, `"median"`, or `"centroid"`. Default is `"complete"`.
#' @param display_numbers Logical. If `TRUE`, numeric values are displayed in the heatmap cells. Default is `FALSE`.
#' @param number_color Character. The color of the numbers displayed in the heatmap cells when `display_numbers = TRUE`. Default is `"black"`.
#' @param color_palette A vector of colors for the heatmap. Default is `c("#05287A", "#DCEEF0")`.
#' @param breaks Numeric. The number of breaks in the color palette. Default is `50`.
#' @param main_title Character. The title of the heatmap. If `NULL`, a default title is generated. Default is `NULL`.
#'
#' @return A distance matrix (as a matrix object) with samples as rows and columns. Row and column names are preserved.
#'
#' @details
#' The function performs the following steps:
#' 1. Checks that the input `count_matrix` is a data frame with row names (sample names).
#' 2. Computes the distance matrix:
#'    - If `use_correlation = TRUE`, the correlation matrix is computed between samples (rows) and converted to a distance matrix using `1 - correlation`.
#'    - If `use_correlation = FALSE`, the distance matrix is computed using the specified distance method.
#' 3. Generates a heatmap of the distance matrix if `plot_heatmap = TRUE`.
#' 4. Returns the distance matrix as a matrix object.
#'
#' @examples
#' # Example count matrix (samples as rows, genes as columns)
#' count_matrix <- data.frame(
#'   Gene1 = c(10, 20, 30),
#'   Gene2 = c(15, 25, 35),
#'   Gene3 = c(12, 22, 32)
#' )
#'
#' # Compute distance matrix using Euclidean distance and plot heatmap
#' dist_matrix <- compute_distance_matrix(count_matrix, use_correlation = FALSE, method = "euclidean")
#'
#' # Compute distance matrix using Spearman correlation and plot heatmap
#' dist_matrix <- compute_distance_matrix(count_matrix, use_correlation = TRUE, correlation_method = "spearman")
#'
#' # Compute distance matrix without plotting the heatmap
#' dist_matrix <- compute_distance_matrix(count_matrix, plot_heatmap = FALSE)
#'
#' # Compute distance matrix and use "average" clustering for the heatmap
#' dist_matrix <- compute_distance_matrix(count_matrix, clustering_method = "average")
#'
#' # Compute distance matrix and display numbers in red
#' dist_matrix <- compute_distance_matrix(count_matrix, display_numbers = TRUE, number_color = "red")
#'
#' @importFrom pheatmap pheatmap
#' @importFrom stats cor dist
#' @export
compute_distance_matrix <- function(count_matrix, use_correlation = FALSE, method = "euclidean", correlation_method = "spearman", plot_heatmap = TRUE,
                                    cluster_rows = TRUE, cluster_cols = TRUE, clustering_method = "complete", display_numbers = FALSE, 
                                    number_color = "black", color_palette = c("#05287A", "#DCEEF0"), breaks = 50, main_title = NULL) {
  # Load necessary libraries
  require(pheatmap)
  require(corrplot)
  
  # Ensure count_matrix is a data frame with samples as rows and genes as columns
  if (!is.data.frame(count_matrix)) {
    stop("Error: count_matrix must be a data frame with samples as row names and genes as columns.")
  }
  
  # Ensure row names are preserved
  if (is.null(rownames(count_matrix))) {
    stop("Error: count_matrix must have row names corresponding to sample names.")
  }
  
  # Compute distance matrix
  if (use_correlation) {
    # Transpose the matrix for correlation (genes as rows, samples as columns)
    cor_matrix <- cor(t(count_matrix), method = correlation_method)
    dist_matrix <- as.dist(1 - cor_matrix)  # Convert correlation to distance
    dist_matrix <- as.matrix(dist_matrix)  # Convert to matrix for visualization
  } else {
    if (method %in% c("euclidean", "manhattan", "canberra", "minkowski")) {
      dist_matrix <- dist(count_matrix, method = method)
      dist_matrix <- as.matrix(dist_matrix)  # Convert to matrix for visualization
    } else {
      stop("Invalid method. Choose from 'euclidean', 'manhattan', 'canberra', or 'minkowski'.")
    }
  }
  
  # Ensure row and column names are preserved
  rownames(dist_matrix) <- rownames(count_matrix)
  colnames(dist_matrix) <- rownames(count_matrix)
  
  # Plot heatmap if required
  if (plot_heatmap) {
    pheatmap(dist_matrix, 
             cluster_rows = cluster_rows,
             cluster_cols = cluster_cols,
             clustering_method = clustering_method,  # Add clustering method
             display_numbers = display_numbers,
             number_color = number_color,  # Add number color
             color = colorRampPalette(color_palette)(breaks),
             main = ifelse(is.null(main_title), paste("Distance Matrix (", ifelse(use_correlation, correlation_method, method), ")"), main_title))
  }
  
  return(dist_matrix)
}
