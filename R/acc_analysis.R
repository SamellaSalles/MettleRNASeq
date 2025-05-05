#' Across-Condition Correlations (ACC) Analysis
#'
#' Calculates correlations between same genes across two different conditions and optionally plots the results.
#' Note: This function only works for comparing exactly 2 conditions at a time.
#'
#' @param data A data frame or matrix containing gene expression data with rows as samples and columns as genes.
#' @param condition_split Numeric vector of length 2 specifying the row indices that separate conditions 
#'        (e.g., c(4, 8) means rows 1-4 are condition 1, rows 5-8 are condition 2).
#' @param method Correlation method to use ("pearson", "kendall", or "spearman"). Default is "spearman".
#' @param plot Logical indicating whether to generate a bar plot. Default is FALSE.
#' @param color Color for the plot elements. Specify either:
#'        - Single color (e.g., "steelblue") for all bars
#'        - Vector of two colors (e.g., c("red", "blue")) for negative/positive correlations
#'        Default is c("red", "steelblue").
#' @param orientation Plot orientation: "vertical" (default) or "horizontal" (uses coord_flip()).
#' @param show_legend Logical indicating whether to show the legend when using two colors. Default is TRUE.
#' @param rotate_labels Logical indicating whether to rotate x-axis labels (vertical plots only). Default is TRUE.
#' @param title Plot title. Default is "Across-Condition Correlations".
#' @param xlab X-axis label. Default is "Gene" (becomes y-label in horizontal orientation).
#' @param ylab Y-axis label. Default is "Correlation" (becomes x-label in horizontal orientation).
#'
#' @return A data frame containing genes and their correlation values, and optionally a ggplot object.
#'
#' @examples
#' \dontrun{
#' # Vertical bars with rotated labels (default)
#' result <- acc_analysis(data, condition_split = c(4, 8), plot = TRUE)
#' 
#' # Horizontal bars (automatically swaps axes)
#' result <- acc_analysis(data, condition_split = c(4, 8), plot = TRUE, 
#'                      orientation = "horizontal")
#' 
#' # Vertical bars without label rotation
#' result <- acc_analysis(data, condition_split = c(4, 8), plot = TRUE,
#'                      rotate_labels = FALSE)
#' }
#' 
#' @importFrom stats cor
#' @importFrom ggplot2 ggplot aes geom_bar theme_minimal labs theme element_text 
#' @importFrom ggplot2 scale_fill_manual coord_flip guide_legend
#' @importFrom dplyr arrange mutate
#' @importFrom rlang .data
#' @export
acc_analysis <- function(data, condition_split, 
                         method = "spearman", 
                         plot = FALSE, 
                         color = c("red", "steelblue"), 
                         orientation = c("vertical", "horizontal"),
                         show_legend = TRUE,
                         rotate_labels = TRUE,
                         title = "Across-Condition Correlations",
                         xlab = "Gene", 
                         ylab = "Correlation") {
  
  # Validate inputs
  if (length(condition_split) != 2) {
    stop("This function only works for comparing exactly 2 conditions. 'condition_split' must be a vector of length 2.")
  }
  
  if (condition_split[2] > nrow(data)) {
    stop("'condition_split' values exceed the number of rows in the data.")
  }
  
  if (length(color) > 2) {
    warning("More than 2 colors provided. Only the first two will be used for negative/positive correlations.")
    color <- color[1:2]
  }
  
  orientation <- match.arg(orientation)
  
  # Check if resamples per condition are equal
  n_cond1 <- condition_split[1]
  n_cond2 <- condition_split[2] - condition_split[1]
  if (n_cond1 != n_cond2) {
    message("Note: Unequal number of resamples per condition (", n_cond1, " vs ", n_cond2, ").",
            "\nFor ideal results, the number of resamples should be the same for both conditions.")
  }
  
  # Convert to data frame if not already
  data <- as.data.frame(data)
  
  # Split data by conditions
  cond1 <- data[1:condition_split[1], , drop = FALSE]
  cond2 <- data[(condition_split[1]+1):condition_split[2], , drop = FALSE]
  
  # Get gene names from column names
  gene_names <- colnames(data)
  if (is.null(gene_names)) {
    gene_names <- paste0("Gene_", 1:ncol(data))
  }
  
  # Calculate correlations for each gene
  correlations <- sapply(1:ncol(data), function(i) {
    cor(cond1[, i], cond2[, i], method = method)
  })
  
  # Create results data frame
  result <- data.frame(
    Gene = factor(gene_names, levels = gene_names),
    Correlation = correlations,
    Direction = ifelse(correlations >= 0, "Positive", "Negative")
  )
  
  # Sort by correlation value
  result <- result %>%
    dplyr::arrange(.data$Correlation) %>%
    dplyr::mutate(Gene = factor(.data$Gene, levels = .data$Gene))
  
  # Generate plot if requested
  if (plot) {
    p <- ggplot2::ggplot(result, ggplot2::aes(x = .data$Gene, y = .data$Correlation)) +
      ggplot2::labs(title = title, x = xlab, y = ylab) +
      ggplot2::theme_minimal()
    
    # Add bars with color handling
    if (length(color) == 1) {
      p <- p + ggplot2::geom_bar(stat = "identity", fill = color)
    } else {
      p <- p + 
        ggplot2::geom_bar(ggplot2::aes(fill = .data$Direction), stat = "identity") +
        ggplot2::scale_fill_manual(values = setNames(color, c("Negative", "Positive")),
                                   guide = if(show_legend) "legend" else "none")
    }
    
    # Apply orientation flip if horizontal
    if (orientation == "horizontal") {
      p <- p + ggplot2::coord_flip()
    } else if (rotate_labels) {
      # Only rotate labels for vertical orientation
      p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    }
    
    print(p)
    return(invisible(list(results = result, plot = p)))
  }
  
  return(result)
}
