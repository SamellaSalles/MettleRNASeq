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
#' @param algorithm Optional algorithm names to use for selection and plotting. With NULL, selection is the union of qualifying resamples across algorithms; plot panels remain separate.
#' @param gene_ids_table Optional gene-ID table from classification. By default, identity metadata on genes_table is used; old tables fall back to positional selection.
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
#'   test_colors = c("#f57600", "green", "red", "purple")
#' )
#' }
#'
#' @export
genes_from_classification_results <- function(eval_table, genes_table, count_matrix,
                                              resamples = NULL, cutoff = NULL,
                                              metric_column = "summary_test_accuracy",
                                              selection_mode = "unique", partition_table = NULL,
                                              plot_type = "summary", plot_resamples_only = TRUE,
                                              return_plot = TRUE, hline_yintercept = 0.7,
                                              hline_color = "black", median_color = "#f57600",
                                              test_colors = NULL, show_legend = TRUE,
                                              algorithm = NULL, gene_ids_table = NULL) {
  selection_mode <- match.arg(selection_mode, c("unique", "intersect"))
  plot_type <- match.arg(plot_type, c("summary", "all", "both"))
  if (!is.data.frame(eval_table) || !is.data.frame(genes_table)) stop("eval_table and genes_table must be data frames.")
  names(eval_table) <- tolower(names(eval_table)); metric_column <- tolower(metric_column)
  if (!"resample" %in% names(eval_table)) stop("eval_table must contain resample.")
  if (!"algorithm" %in% names(eval_table)) eval_table$algorithm <- "classification"
  if (!is.null(algorithm)) {
    if (!all(algorithm %in% eval_table$algorithm)) stop("Unknown algorithm in eval_table.")
    eval_table <- eval_table[eval_table$algorithm %in% algorithm, , drop = FALSE]
  }
  if (is.null(resamples)) {
    if (is.null(cutoff) || length(cutoff) != 1L || !is.numeric(cutoff) || !is.finite(cutoff)) stop("Supply resamples or a finite cutoff.")
    if (!metric_column %in% names(eval_table)) stop("Metric column not found in eval_table.")
    keep <- is.finite(eval_table[[metric_column]]) & eval_table[[metric_column]] >= cutoff
    if ("status" %in% names(eval_table)) keep <- keep & eval_table$status == "ok"
    selected_resamples <- unique(eval_table$resample[keep])
  } else {
    if (!is.numeric(resamples) || anyNA(resamples) || any(!is.finite(resamples)) ||
        any(resamples < 1L | resamples != floor(resamples))) stop("resamples must contain positive integer IDs.")
    selected_resamples <- unique(resamples)
  }
  if (!length(selected_resamples)) stop("No resamples passed the condition.")
  columns <- paste0("Resample_", selected_resamples)
  if (!all(columns %in% names(genes_table))) stop("Some selected resamples not found in genes_table.")
  x <- .mettle_numeric_matrix(count_matrix, "count_matrix")
  if (is.null(gene_ids_table)) gene_ids_table <- attr(genes_table, "gene_ids")
  combine <- function(values) if (selection_mode == "unique") unique(unlist(values, use.names = FALSE)) else Reduce(intersect, values)
  if (!is.null(gene_ids_table)) {
    if (!is.data.frame(gene_ids_table) || !all(columns %in% names(gene_ids_table))) stop("gene_ids_table must contain every selected resample.")
    ids <- as.character(combine(lapply(gene_ids_table[, columns, drop = FALSE], as.character)))
    .mettle_names(colnames(x), "Gene column names")
    if (anyNA(ids) || !all(ids %in% colnames(x))) stop("Selected gene IDs do not match count_matrix columns.")
    positions <- match(ids, colnames(x))
  } else {
    # Compatibility for old positional tables, which have no identity metadata.
    positions <- combine(as.list(genes_table[, columns, drop = FALSE]))
    if (!is.numeric(positions) || anyNA(positions) || any(!is.finite(positions)) ||
        any(positions < 1 | positions != floor(positions) | positions > ncol(x))) stop("Invalid gene positions for count_matrix.")
    ids <- colnames(x)[positions]
  }
  filtered <- as.data.frame(x[, positions, drop = FALSE])
  plot_obj <- NULL
  if (return_plot) {
    plot_data <- if (plot_resamples_only) eval_table[eval_table$resample %in% selected_resamples, , drop = FALSE] else eval_table
    if (plot_type %in% c("all", "both")) {
      if (is.null(partition_table) || !all(c("Train_Splits", "Test_Splits") %in% names(partition_table))) stop("partition_table is required for individual test accuracies.")
      if (!"data_split" %in% names(plot_data) || anyNA(plot_data$data_split) ||
          any(plot_data$data_split < 1 | plot_data$data_split > nrow(partition_table))) stop("Invalid data_split values.")
      plot_data$Partition_Label <- paste0("Train: [", partition_table$Train_Splits[plot_data$data_split],
        "] | Test: [", partition_table$Test_Splits[plot_data$data_split], "]")
    }
    p <- ggplot2::ggplot()
    if (plot_type %in% c("all", "both")) {
      if (!"test_accuracy" %in% names(plot_data)) stop("test_accuracy is required for individual plots.")
      p <- p + ggplot2::geom_point(data = plot_data,
        ggplot2::aes(x = factor(.data$resample), y = .data$test_accuracy, color = .data$Partition_Label), size = 2, alpha = .7)
    }
    if (plot_type %in% c("summary", "both")) {
      if (!"summary_test_accuracy" %in% names(plot_data)) stop("summary_test_accuracy is required for summary plots.")
      summary_df <- unique(plot_data[, c("algorithm", "resample", "summary_test_accuracy"), drop = FALSE])
      if (anyDuplicated(summary_df[,c("algorithm", "resample")])) stop("Inconsistent summary values for an algorithm/resample pair.")
      summary_df <- summary_df[order(summary_df$algorithm, summary_df$resample), , drop = FALSE]
      p <- p + ggplot2::geom_line(data = summary_df,
        ggplot2::aes(x = factor(.data$resample), y = .data$summary_test_accuracy, group = .data$algorithm), color = median_color, linewidth = 1) +
        ggplot2::geom_point(data = summary_df,
          ggplot2::aes(x = factor(.data$resample), y = .data$summary_test_accuracy), color = median_color, size = 2)
    }
    p <- p + ggplot2::geom_hline(yintercept = hline_yintercept, linetype = "dashed", color = hline_color, linewidth = .6) +
      ggplot2::labs(x = "Resample", y = "Accuracy", title = "Classification Accuracy by Resample") + ggplot2::theme_bw()
    if (length(unique(plot_data$algorithm)) > 1L) p <- p + ggplot2::facet_wrap(~algorithm)
    if (plot_type %in% c("all", "both")) {
      labels <- unique(plot_data$Partition_Label)
      colors <- if (is.null(test_colors)) grDevices::hcl.colors(length(labels), "Dark 3") else test_colors
      if (is.null(names(colors))) {
        if (length(colors) < length(labels)) stop("Not enough test_colors for partitions.")
        colors <- colors[seq_along(labels)]
        names(colors) <- labels
      }
      p <- p + ggplot2::scale_color_manual(name = "Partition", values = colors)
    }
    if (!show_legend) p <- p + ggplot2::theme(legend.position = "none")
    print(p); plot_obj <- p
  }
  list(selected_resamples = selected_resamples, selected_gene_positions = positions,
    filtered_matrix = filtered, plot = plot_obj, selected_gene_ids = ids)
}
