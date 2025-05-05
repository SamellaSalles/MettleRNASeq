#' Intra-Condition Correlation (ICC) Analysis
#'
#' Calculates gene-gene correlations within each condition. Returns proper correlation matrices.
#' Note: Genes must be in columns (variables) and samples in rows (observations).
#'
#' @param data A data frame or matrix with samples in rows and genes in columns.
#' @param condition_split Numeric vector of length 2 specifying the row indices that separate conditions.
#' @param method Correlation method ("pearson", "kendall", or "spearman"). Default is "spearman".
#' @param plot Logical indicating whether to generate correlation heatmaps. Default is FALSE.
#' @param colors Vector of colors for the heatmap gradient. Default is c("#f57600", "white", "#0073e6").
#' @param plot_type Type of heatmap ("full", "upper", or "lower"). Default is "lower".
#' @param show_coef Logical indicating whether to show correlation coefficients. Default is FALSE.
#' @param tl.cex Text size for gene names. Default is 0.7.
#' @param number.cex Text size for correlation coefficients. Default is 0.7.
#' @param order Method to order the correlation matrix ("AOE", "FPC", "hclust", "alphabet"). Default is "AOE".
#' 
#' @return A list containing proper matrix objects:
#' \itemize{
#'   \item correlation_matrix_cond1 - Gene x gene correlation matrix for condition 1 (class "matrix")
#'   \item correlation_matrix_cond2 - Gene x gene correlation matrix for condition 2 (class "matrix")
#'   \item plots - (Optional) List of heatmap plots when plot=TRUE
#' }
#'
#' @examples
#' \dontrun{
#' data <- matrix(rnorm(80), nrow=8, ncol=10)  # 8 samples (4+4), 10 genes
#' colnames(data) <- paste0("Gene_", 1:10)
#' icc_result <- icc_analysis(data, condition_split = c(4, 8))
#' class(icc_result$correlation_matrix_cond1)  # "matrix" "array" 
#' }
#' 
#' @importFrom stats cor
#' @importFrom corrplot corrplot
#' @importFrom grDevices colorRampPalette
#' @export
icc_analysis <- function(data, condition_split, 
                         method = "spearman",
                         plot = FALSE,
                         colors = c("#f57600", "white", "#0073e6"),
                         plot_type = c("lower", "upper", "full"),
                         show_coef = FALSE,
                         tl.cex = 0.7,
                         number.cex = 0.7,
                         order = c("AOE", "FPC", "hclust", "alphabet")) {
  
  # Validate input structure
  if (nrow(data) <= 1 || ncol(data) <= 1) {
    stop("Data must have >1 row (samples) and >1 column (genes). Genes must be in columns.")
  }
  
  if (length(condition_split) != 2) {
    stop("condition_split must be a vector of length 2 specifying the row indices separating conditions.")
  }
  
  if (condition_split[2] != nrow(data)) {
    stop("Second value in condition_split must match total number of samples (nrow(data)).")
  }
  
  plot_type <- match.arg(plot_type)
  order <- match.arg(order)
  
  # Convert to matrix (preserving dimnames)
  data <- as.matrix(data)
  gene_names <- colnames(data)
  
  if (is.null(gene_names)) {
    gene_names <- paste0("Gene_", 1:ncol(data))
    colnames(data) <- gene_names
  }
  
  # Split data by conditions
  cond1 <- data[1:condition_split[1], , drop = FALSE]
  cond2 <- data[(condition_split[1]+1):condition_split[2], , drop = FALSE]
  
  # Calculate correlation matrices (genes x genes) - ensures matrix output
  corr_cond1 <- matrix(cor(cond1, method = method), 
                       ncol = ncol(cond1),
                       dimnames = list(gene_names, gene_names))
  
  corr_cond2 <- matrix(cor(cond2, method = method), 
                       ncol = ncol(cond2),
                       dimnames = list(gene_names, gene_names))
  
  # Generate plots if requested
  plot_list <- list()
  if (plot) {
    color_palette <- colorRampPalette(colors)(100)
    
    plot_list$plot_cond1 <- corrplot::corrplot(
      corr_cond1,
      method = "color",
      col = color_palette,
      type = plot_type,
      tl.cex = tl.cex,
      tl.col = "black",
      number.cex = number.cex,
      order = order,
      addCoef.col = if(show_coef) "black" else NULL,
      mar = c(0, 0, 2, 0),
      title = "Condition 1"
    )
    
    plot_list$plot_cond2 <- corrplot::corrplot(
      corr_cond2,
      method = "color",
      col = color_palette,
      type = plot_type,
      tl.cex = tl.cex,
      tl.col = "black",
      number.cex = number.cex,
      order = order,
      addCoef.col = if(show_coef) "black" else NULL,
      mar = c(0, 0, 2, 0),
      title = "Condition 2"
    )
  }
  
  # Prepare results - ensure we're returning proper matrices
  results <- list(
    correlation_matrix_cond1 = corr_cond1,
    correlation_matrix_cond2 = corr_cond2
  )
  
  if (plot) {
    results$plots <- plot_list
  }
  
  return(results)
}
