#' Advanced Classification with Gene Resampling and Replicate Splitting
#'
#' This function performs classification using selected algorithms across multiple resamples of genes,
#' with all possible train/test replicate partitions. It outputs evaluation metrics, selected genes,
#' partitioning details, and trained models.
#'
#' @param input_data A data frame output from \code{prepare_input_matrix()}, with genes as rows,
#' replicates (e.g., \code{Rep_1}, \code{Rep_2}, ...) as columns, and a \code{Condition} column.
#' @param algorithms A character vector of classification algorithms (e.g., \code{"rf"}, \code{"ranger"}, \code{"svmLinear2"}).
#' @param hyper_params A named list of tuning grids (data.frames) for each algorithm. Each list element must contain
#' the hyperparameters specific to the corresponding algorithm. See details below.
#' @param num_resamples Number of bootstrap resamples (i.e., how many different sets of genes to sample).
#' @param num_genes Number of genes to randomly sample from the first condition in each resample.
#' @param train_percentage Percentage of replicates to use for training (e.g., 0.75).
#' @param metrics A character vector of evaluation metrics (currently only \code{"Accuracy"} is supported).
#' @param summary_stat Statistic to summarize test performance per resample. Either \code{"median"} or \code{"mean"}.
#' @param replace_genes Logical; whether to sample genes with replacement.
#' @param parallel Logical; whether to allow parallel execution during model training (if supported by the algorithm).
#' @param seed Integer; random seed for reproducibility.
#'
#' @details
#' Genes are sampled only from the first condition, then matched across all conditions. For each resample, 
#' the function uses all possible combinations of train/test replicate splits based on the given training percentage.
#'
#' The \code{hyper_params} argument must be a named list where each name corresponds to an algorithm (e.g., \code{"rf"} or \code{"ranger"}), and each value is a data frame with hyperparameter settings. Example:
#'
#' \preformatted{
#' hyper_params <- list(
#'   rf = data.frame(mtry = c(1, 5, 10, sqrt(length(cols)))), 
#'   ranger = data.frame(mtry = c(1, 5, 10, sqrt(length(cols))), 
#'                       splitrule = "extratrees", 
#'                       min.node.size = c(2, 3))
#' )
#' }
#'
#' Output metrics include training and test accuracy per split and per resample, with summaries per resample.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{eval_table}{A data frame with algorithm name, resample number, data split, train/test replicates, and metrics.}
#'   \item{partition_table}{A data frame listing all unique train/test replicate combinations used.}
#'   \item{genes_table}{A data frame where each column is a resample and each row is the position (index) of a gene selected from the first condition.}
#'   \item{models}{A nested list of trained models per algorithm, resample, and data split.}
#' }
#'
#' @examples
#' \dontrun{
#' input_data <- prepare_input_matrix(my_counts, sample_vector = sample_names, condition_vector = condition_labels)
#' hyper_params <- list(
#'   rf = data.frame(mtry = c(1, 5, 10))
#' )
#' results <- classification_approach_advanced(
#'   input_data = input_data,
#'   algorithms = c("rf"),
#'   hyper_params = hyper_params,
#'   num_resamples = 5,
#'   num_genes = 50
#' )
#' }
#'
#' @export
classification_approach_advanced <- function(input_data, algorithms, hyper_params, num_resamples, num_genes,
                                             train_percentage = 0.75, metrics = "Accuracy", summary_stat = "median",
                                             replace_genes = TRUE, parallel=FALSE, seed=123) {
  
  library(caret)
  library(dplyr)
  library(e1071)
  library(pROC)
  
  # Ensure input columns
  if (!all(c("Gene", "Condition") %in% colnames(input_data))) {
    stop("Input data must have 'Gene' and 'Condition' columns.")
  }
  
  # Detect replicate columns dynamically
  replicate_cols <- setdiff(colnames(input_data), c("Gene", "Condition"))
  num_replicates <- length(replicate_cols)
  num_conditions <- length(unique(input_data$Condition))
  genes_per_condition <- nrow(input_data) / num_conditions
  
  # Generate train/test splits (fixed partitions)
  train_test_splits <- combn(1:num_replicates, round(train_percentage * num_replicates), simplify = FALSE)
  train_test_splits <- lapply(train_test_splits, function(train_indices) {
    list(
      train_replicates = train_indices,
      test_replicates = setdiff(1:num_replicates, train_indices)
    )
  })
  
  # Initialize results
  eval_table <- data.frame()
  partition_table <- data.frame()
  genes_table <- list()  # Now we store genes in a list for each resample
  models <- list()
  
  # Create resamples
  set.seed(seed)
  resample_indices <- lapply(1:num_resamples, function(b) {
    selected_genes <- sample(1:genes_per_condition, num_genes, replace = replace_genes)
    gene_indices <- unlist(lapply(0:(num_conditions - 1), function(i) {
      selected_genes + i * genes_per_condition
    }))
    list(gene_indices = gene_indices, train_test_splits = train_test_splits)
  })
  names(resample_indices) <- paste0("Bootstrap", 1:num_resamples)
  
  # Build genes_table as a clean data frame
  genes_table <- as.data.frame(do.call(cbind, lapply(resample_indices, function(x) x$gene_indices)))[1:num_genes,]
  colnames(genes_table) <- paste0("Resample_", seq_along(resample_indices))
  
  # Create a partition table (fixed partitions)
  partition_counter <- 1
  for (split_idx in seq_along(train_test_splits)) {
    split_info <- train_test_splits[[split_idx]]
    partition_table <- rbind(partition_table, data.frame(
      Partition_ID = paste("Partition", partition_counter),
      Train_Splits = paste(replicate_cols[split_info$train_replicates], collapse = ","),
      Test_Splits = paste(replicate_cols[split_info$test_replicates], collapse = ",")
    ))
    partition_counter <- partition_counter + 1
  }
  
  # Iterate over classification algorithms
  for (m in algorithms) {
    accuracy_df <- data.frame()
    models_boot <- list()
    
    # Iterate over resamples
    for (i in seq_along(resample_indices)) {
      resample_info <- resample_indices[[i]]
      gene_subset <- input_data[resample_info$gene_indices, ]
      cols_df <- data.frame()
      
      # Iterate over train/test splits
      for (split_idx in seq_along(resample_info$train_test_splits)) {
        split_info <- resample_info$train_test_splits[[split_idx]]
        train_cols_idx <- split_info$train_replicates+1
        train_data <- gene_subset[, c(1,train_cols_idx,6)]
        indices_train_data <- 1:nrow(train_data)
        train_data$Condition <- factor(train_data$Condition)
        train_cols <- replicate_cols[split_info$train_replicates]
        
        test_cols_idx <- split_info$test_replicates+1
        test_data <- gene_subset[,!(1:length(gene_subset)) %in% c(train_cols_idx,6)]
        test_data <- cbind(test_data,replicate(length(train_cols_idx)-1,test_data[2]))
        test_data$Condition <- gene_subset$Condition
        colnames(test_data) <- c("Gene",names(gene_subset[train_cols_idx]),"Condition")
        test_cols <- replicate_cols[split_info$test_replicates]
        
        # Set up trainControl for bootstrapping and resampling
        trControl <- trainControl(
          method = "boot", 
          number = 1,
          index = list(indices_train_data),
          indexOut = list(indices_train_data)
        )
        
        # Define the tuning grid (single row for non-resampling or multi rows if tuning)
        tuneGrid <- hyper_params[[m]]
        
        set.seed(seed)
        result <- tryCatch({
          model <- train(
            Condition~.,
            data = train_data,
            method = m,
            metric = "Accuracy",
            preProcess = "center",
            trControl = trControl,
            tuneGrid = tuneGrid,
            allowParallel= parallel
          )
          
          acc_treino <- model$results%>%
            filter(mtry==model$bestTune$mtry)%>%
            select(Accuracy)
          
          predictions <- predict(model, newdata = test_data[1:5])
          ref <- factor(test_data$Condition)
          metrics_list <- list()
          
          conf_matrix <- confusionMatrix(predictions, ref)
          metrics_list$accuracy <- conf_matrix$overall["Accuracy"]
          #metrics_list$accuracy <- mean(predictions == ref)
          
          list(model = model, metrics = metrics_list)
        }, error = function(e) return(NULL))
        
        if (is.null(result) || length(result$metrics) < 1) next
        
        # Collect metrics for the current resample and split
        metrics_row <- data.frame(
          algorithm = m,
          resample = i,
          data_split = split_idx,
          train_partition = paste(train_cols, collapse = ","),
          test_partition = paste(test_cols, collapse = ",")
        )
        
        for (metric in names(result$metrics)) {
          metrics_row[[paste0("train_accuracy")]] <- acc_treino$Accuracy
          metrics_row[[paste0("test_", metric)]] <- result$metrics[[metric]]
        }
        
        # Add to eval table
        cols_df <- rbind(cols_df, metrics_row)
        models_boot[[paste0("resample_", i, "_split_", split_idx)]] <- result$model
      }
      
      if (nrow(cols_df) > 0) {
        for (metric in names(result$metrics)) {
          cols_df[[paste0("summary_test_", metric)]] <- if (summary_stat == "median") {
            median(cols_df[[paste0("test_", metric)]], na.rm = TRUE)
          } else {
            mean(cols_df[[paste0("test_", metric)]], na.rm = TRUE)
          }
        }
      }
      
      accuracy_df <- rbind(accuracy_df, cols_df)
      
    }
    
    eval_table <- rbind(eval_table, accuracy_df)
    models[[m]] <- models_boot
  }
  
  # Return genes_table as a list of gene positions for each resample
  return(list(
    eval_table = eval_table,
    partition_table = partition_table,  # return partition table with Partition_ID
    genes_table = genes_table,  # List of gene positions for each resample
    models = models
  ))
}

