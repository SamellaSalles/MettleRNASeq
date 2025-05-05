#' Absolute Difference Between Intra-Condition Correlations (ICC)
#'
#' Computes absolute differences between gene-gene correlation matrices from two conditions,
#' returns results in a tidy format, and provides visualization options.
#'
#' @param icc_result List containing ICC results from icc_analysis()
#' @param cutoff Optional numeric threshold for filtering differences (NULL returns all)
#' @param plot Logical indicating whether to generate a heatmap (default: FALSE)
#' @param heatmap_colors Colors for heatmap gradient (default: yellow to red)
#' @param heatmap_breaks Break points for color bins (default: seq(0, 2, by=0.5))
#' @param show_rownames Logical indicating whether to show row names in heatmap (default: TRUE)
#' @param show_colnames Logical indicating whether to show column names in heatmap (default: TRUE)
#' @param cluster_rows Logical indicating whether to cluster rows (default: TRUE)
#' @param cluster_cols Logical indicating whether to cluster columns (default: TRUE)
#' @param fontsize Base fontsize for heatmap (default: 8)
#'
#' @return A list containing:
#' \itemize{
#'   \item absolute_diff_matrix - Full matrix of absolute differences
#'   \item difference_table - Data frame of gene pairs and their differences
#'   \item heatmap - (Optional) Heatmap object when plot=TRUE
#' }
#'
#' @examples
#' \dontrun{
#' data <- matrix(rnorm(80), nrow=8, ncol=10)
#' colnames(data) <- paste0("Gene_", 1:10)
#' icc_result <- icc_analysis(data, condition_split = c(4, 8))
#'
#' # Basic usage with gene names visible by default
#' diff_result <- icc_difference(icc_result, plot = TRUE)
#' }
#' 
#' @importFrom pheatmap pheatmap
#' @export
icc_difference <- function(icc_result, 
                           cutoff = NULL,
                           plot = FALSE,
                           heatmap_colors = colorRampPalette(c("#fffd8d", "yellow", "red")),
                           heatmap_breaks = seq(0, 2, by = 0.5),
                           show_rownames = TRUE,
                           show_colnames = TRUE,
                           cluster_rows = TRUE,
                           cluster_cols = TRUE,
                           fontsize = 8) {
  
  # Validate input
  if (!all(c("correlation_matrix_cond1", "correlation_matrix_cond2") %in% names(icc_result))) {
    stop("Input must be result from icc_analysis() containing correlation matrices for both conditions")
  }
  
  corr1 <- icc_result$correlation_matrix_cond1
  corr2 <- icc_result$correlation_matrix_cond2
  
  if (!identical(dim(corr1), dim(corr2)) || !identical(rownames(corr1), rownames(corr2))) {
    stop("Correlation matrices must have identical dimensions and gene names")
  }
  
  # Compute absolute differences
  abs_diff <- abs(corr1 - corr2)
  
  # Apply cutoff if specified
  if (!is.null(cutoff)) {
    filtered_diff <- abs_diff
    filtered_diff[filtered_diff < cutoff] <- NA
    indices <- which(!is.na(filtered_diff), arr.ind = TRUE)
  } else {
    indices <- which(!is.na(abs_diff), arr.ind = TRUE)
  }
  
  # Create difference table
  diff_table <- data.frame(
    Gene1 = rownames(abs_diff)[indices[, 1]],
    Gene2 = colnames(abs_diff)[indices[, 2]],
    Corr_Condition1 = corr1[indices],
    Corr_Condition2 = corr2[indices],
    Abs_Difference = abs_diff[indices],
    stringsAsFactors = FALSE
  )
  
  # Sort by absolute difference (descending)
  diff_table <- diff_table[order(-diff_table$Abs_Difference), ]
  rownames(diff_table) <- NULL
  
  # Generate heatmap if requested
  heatmap_obj <- NULL
  if (plot) {
    plot_mat <- abs_diff
    if (!is.null(cutoff)) {
      plot_mat[abs_diff < cutoff] <- NA
    }
    
    heatmap_obj <- pheatmap::pheatmap(
      mat = plot_mat,
      color = heatmap_colors(length(heatmap_breaks) - 1),
      breaks = heatmap_breaks,
      na_col = "white",
      show_rownames = show_rownames,
      show_colnames = show_colnames,
      cluster_rows = cluster_rows,
      cluster_cols = cluster_cols,
      fontsize = fontsize,
      silent = TRUE
    )
    
    grid::grid.newpage()
    grid::grid.draw(heatmap_obj$gtable)
  }
  
  # Return results
  result <- list(
    absolute_diff_matrix = abs_diff,
    difference_table = diff_table
  )
  
  if (plot) {
    result$heatmap <- heatmap_obj
  }
  
  return(result)
}
