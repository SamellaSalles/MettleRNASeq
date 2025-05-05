#' Prepare Input Matrix for Classification
#'
#' This function prepares an input matrix for classification by organizing count data
#' into a wide format with consistent replicate naming.
#'
#' @param count_matrix A matrix or data frame of count data, where rows are genes and columns are samples.
#' @param design_table A data frame containing sample and condition information. Default is NULL.
#' @param sample_col The name of the column in `design_table` containing sample names. Default is NULL.
#' @param condition_col The name of the column in `design_table` containing condition names. Default is NULL.
#' @param sample_vector A vector of sample names. Default is NULL.
#' @param condition_vector A vector of condition names corresponding to `sample_vector`. Default is NULL.
#'
#' @return A data frame in wide format with genes as rows, replicates as columns (named "Rep1", "Rep2", etc.),
#'         and conditions as the last column.
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
prepare_input_matrix <- function(count_matrix, 
                                 design_table = NULL, 
                                 sample_col = NULL, 
                                 condition_col = NULL,
                                 sample_vector = NULL, 
                                 condition_vector = NULL) {
  
  # Validate input matrix
  if (ncol(count_matrix) > nrow(count_matrix)) {
    count_matrix <- as.data.frame(t(count_matrix))
  }
  count_matrix <- as.data.frame(count_matrix)
  
  # Create condition mapping
  if (!is.null(design_table)) {
    if (is.null(sample_col) || is.null(condition_col)) {
      stop("When using design_table, both sample_col and condition_col must be specified.")
    }
    condition_map <- design_table %>%
      dplyr::select(Sample = !!sym(sample_col), 
                    Condition = !!sym(condition_col))
  } else if (!is.null(sample_vector) && !is.null(condition_vector)) {
    if (length(sample_vector) != length(condition_vector)) {
      stop("sample_vector and condition_vector must have equal length.")
    }
    condition_map <- data.frame(Sample = sample_vector, 
                                Condition = condition_vector, 
                                stringsAsFactors = FALSE)
  } else {
    stop("Must provide either design_table or both sample_vector and condition_vector")
  }
  
  # Check sample matching
  if (!all(condition_map$Sample %in% colnames(count_matrix))) {
    stop("All samples in condition mapping must exist in count_matrix columns.")
  }
  
  # Verify balanced replicates
  rep_counts <- table(condition_map$Condition)
  if (length(unique(rep_counts)) > 1) {
    stop("All conditions must have equal replicates. Found: ", 
         paste(rep_counts, collapse = ", "))
  }
  
  # Process condition ordering
  if (!is.null(condition_vector)) {
    condition_map$Condition <- factor(condition_map$Condition, 
                                      levels = unique(condition_vector))
  }
  
  # Prepare long format data
  count_matrix$Gene <- rownames(count_matrix)
  
  long_data <- count_matrix %>%
    tidyr::pivot_longer(cols = -Gene, 
                        names_to = "Sample", 
                        values_to = "Expression") %>%
    dplyr::left_join(condition_map, by = "Sample") %>%
    dplyr::arrange(Condition, match(Gene, rownames(count_matrix))) %>%
    dplyr::group_by(Condition, Gene) %>%
    dplyr::mutate(Replicate = paste0("Rep_", row_number())) %>%
    dplyr::ungroup() %>%
    dplyr::select(Gene, Replicate, Expression, Condition)
  
  # Convert to wide format with clean column names
  wide_data <- long_data %>%
    tidyr::pivot_wider(names_from = Replicate, 
                       values_from = Expression) %>%
    dplyr::relocate(Condition, .after = last_col())
  
  return(wide_data)
}

