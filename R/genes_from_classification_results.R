#' Select Genes Based on Classification Results and Visualize Accuracy
#'
#' This function filters genes from a count matrix based on resamples selected from classification results
#' (`eval_table`) using a specified accuracy metric or vector of resample IDs. It also provides a customizable
#' accuracy plot showing test and/or summary performance across resamples.
#'
#' @param eval_table A data frame returned from the classification function, containing evaluation metrics
#' for each resample and data split.
#' @param genes_table A data frame (output of the classification function) containing the gene positions
#' used in each resample, named like `"Resample_1"`, `"Resample_2"`, etc.
#' @param count_matrix A count matrix with genes in the columns. This will be filtered to keep only
#' the selected gene columns.
#' @param resamples Optional vector of resample IDs to select (e.g., \code{c(1, 3, 5)}). If NULL, \code{cutoff}
#' will be used.
#' @param cutoff Optional numeric value used to filter resamples based on the specified \code{metric_column}.
#' Only one of \code{resamples} or \code{cutoff} is required.
#' @param metric_column Character. The name of the metric column to use for filtering resamples when \code{cutoff}
#' is provided. Typically "summary_test_accuracy" or "test_accuracy".
#' @param selection_mode Either \code{"unique"} to take the union of selected genes across resamples,
#' or \code{"intersect"} to take only the genes common to all selected resamples.
#' @param partition_table Optional data frame returned by the classification function, showing
#' which replicates were used in each data split (used for coloring test accuracy points).
#' @param plot_type Type of accuracy plot to generate. Options are:
#' \code{"summary"} (median/mean test accuracy),
#' \code{"all"} (individual test accuracies),
#' or \code{"both"} (both overlays).
#' @param plot_resamples_only Logical. If TRUE (default), only the selected resamples are plotted.
#' If FALSE, all resamples from \code{eval_table} are shown.
#' @param return_plot Logical. If TRUE (default), displays and returns the ggplot object.
#' @param hline_yintercept Numeric value for drawing a horizontal threshold line (e.g., accuracy = 0.7).
#' @param hline_color Color of the horizontal reference line (default = "black").
#' @param median_color Color used for the median/mean accuracy line and points (default = "#f57600").
#' @param test_colors Optional named vector of colors to manually assign colors to test accuracy points,
#' based on the partitions (Train/Test splits).
#' @param show_legend Logical. Whether to show the plot legend (default = TRUE).
#'
#' @return A list with:
#' \itemize{
#'   \item \code{selected_resamples}: The resamples used to select genes.
#'   \item \code{selected_gene_positions}: The positions of the selected genes (column indices of \code{count_matrix}).
#'   \item \code{filtered_matrix}: A filtered count matrix (as a \code{data.frame}) containing only the selected genes.
#'   \item \code{plot}: The \code{ggplot2} object generated (or NULL if \code{return_plot = FALSE}).
#' }
#'
#' @examples
#' \dontrun{
#' result <- genes_from_classification_results(
#'   eval_table = classification$eval_table,
#'   genes_table = classification$genes_table,
#'   count_matrix = my_counts,
#'   partition_table = classification$partition_table,
#'   cutoff = 0.6,
#'   metric_column = "summary_test_accuracy",
#'   selection_mode = "unique",
#'   plot_type = "both",
#'   plot_resamples_only = TRUE,
#'   hline_yintercept = 0.7,
#'   median_color = "#f57600",
#'   test_colors = c("#f57600", "green", "red", "purple"))
#' )
#' }
#'
#' @export
genes_from_classification_results <- function(eval_table,
                                              genes_table,
                                              count_matrix,
                                              resamples = NULL,
                                              cutoff = NULL,
                                              metric_column = "summary_test_accuracy",
                                              selection_mode = "unique",
                                              partition_table = NULL,
                                              plot_type = "summary",  # "all", "summary", or "both"
                                              plot_resamples_only = TRUE,
                                              return_plot = TRUE,
                                              hline_yintercept = 0.7,
                                              hline_color = "black",
                                              median_color = "#f57600",
                                              test_colors = NULL,
                                              show_legend = TRUE) {
  library(dplyr)
  library(ggplot2)
  
  # --- Validations ---
  if (is.null(resamples) && is.null(cutoff)) {
    stop("Provide either resamples or a cutoff.")
  }
  
  colnames(eval_table) <- tolower(colnames(eval_table))
  metric_column <- tolower(metric_column)
  
  if (!(metric_column %in% colnames(eval_table))) {
    stop(paste("Column", metric_column, "not found in eval_table."))
  }
  
  if (plot_type %in% c("all", "both") && is.null(partition_table)) {
    stop("partition_table is required when plotting individual test accuracies.")
  }
  
  # --- Get resamples ---
  selected_resamples <- if (!is.null(resamples)) {
    unique(resamples)
  } else {
    unique(eval_table$resample[eval_table[[metric_column]] >= cutoff])
  }
  
  if (length(selected_resamples) == 0) {
    stop("No resamples passed the condition.")
  }
  
  # --- Extract gene positions ---
  selected_columns <- paste0("Resample_", selected_resamples)
  if (!all(selected_columns %in% colnames(genes_table))) {
    stop("Some selected resamples not found in genes_table.")
  }
  
  selected_gene_matrix <- genes_table[, selected_columns, drop = FALSE]
  
  selected_positions <- if (selection_mode == "unique") {
    unique(unlist(selected_gene_matrix))
  } else if (selection_mode == "intersect") {
    Reduce(intersect, as.list(selected_gene_matrix))
  } else {
    stop("selection_mode must be 'unique' or 'intersect'")
  }
  
  # --- Filter count_matrix (columns only) ---
  if (any(selected_positions > ncol(count_matrix))) {
    stop("Some gene positions exceed column count of count_matrix.")
  }
  
  filtered_matrix <- as.matrix(count_matrix[, selected_positions, drop = FALSE])
  
  # --- Plotting ---
  plot_obj <- NULL
  if (return_plot) {
    plot_data <- if (plot_resamples_only) {
      eval_table %>% filter(resample %in% selected_resamples)
    } else {
      eval_table
    }    
    # Join with partition_table to get test split IDs
    if (!is.null(partition_table)) {
      partition_table <- partition_table %>%
        mutate(data_split = row_number(),
               Partition_Label = paste0("Train: [", Train_Splits, "] | Test: [", Test_Splits, "]"))
      plot_data <- left_join(plot_data, partition_table, by = "data_split")
    }
    
    # Build plot
    p <- ggplot()
    
    # --- Plot individual test points by partition split ---
    if (plot_type %in% c("all", "both") && "test_accuracy" %in% colnames(plot_data)) {
      p <- p + geom_point(aes(x = factor(resample),
                              y = test_accuracy,
                              color = Partition_Label),
                          data = plot_data,
                          size = 2,
                          alpha = 0.7)
    }
    
    # --- Plot summary accuracy per resample ---
    if (plot_type %in% c("summary", "both") && "summary_test_accuracy" %in% colnames(plot_data)) {
      summary_df <- plot_data %>%
        group_by(resample) %>%
        summarise(summary = median(test_accuracy, na.rm = TRUE))
      
      p <- p +
        geom_line(data = summary_df, aes(x = factor(resample), y = summary, group = 1),
                  color = median_color, linewidth = 1) +
        geom_point(data = summary_df, aes(x = factor(resample), y = summary),
                   color = median_color, size = 2)
    }
    
    p <- p +
      geom_hline(yintercept = hline_yintercept,
                 linetype = "dashed", color = hline_color, size = 0.6) +
      labs(x = "Resample", y = "Accuracy", title = "Classification Accuracy by Resample") +
      theme_bw()
    
    if (plot_type %in% c("all", "both") && !is.null(partition_table)) {
      if (!is.null(test_colors)) {
        p <- p + scale_color_manual(name = "Partition", values = test_colors)
      } else {
        p <- p + scale_color_brewer(name = "Partition", palette = "Set1")
      }
    }
    
    if (!show_legend) {
      p <- p + theme(legend.position = "none")
    }
    
    print(p)
    plot_obj <- p
  }
  
  # --- Return ---
  return(list(
    selected_resamples = selected_resamples,
    selected_gene_positions = selected_positions,
    filtered_matrix = as.data.frame(filtered_matrix),
    plot = plot_obj
  ))
}
