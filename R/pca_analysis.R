#' Perform PCA Analysis and Generate Plots
#'
#' This function performs Principal Component Analysis (PCA) on a count matrix and provides options for visualizing the results.
#' It computes eigenvalues, variance, cumulative variance, and cos2 values to assess how well samples are represented by each principal component.
#' The function also generates 2D and 3D PCA plots, with optional coloring and shaping based on specified columns of a design matrix or direct vectors.
#'
#' @param count_matrix A data frame or matrix where rows are samples and columns are genes. Row names should be sample names.
#' @param design_matrix A data frame specifying the conditions, batches, or other metadata for each sample. Default is `NULL`.
#' @param color_by Character or vector. If a character, it should be the column name in `design_matrix` to use for coloring samples. If a vector, it should be a factor or numeric vector for coloring. Default is `NULL`.
#' @param shape_by Character or vector. If a character, it should be the column name in `design_matrix` to use for shaping samples. If a vector, it should be a factor or numeric vector for shaping. Default is `NULL`.
#' @param center Logical. Whether to center the data before performing PCA. Default is `TRUE`.
#' @param scale Logical. Whether to scale the data before performing PCA. Default is `FALSE`.
#' @param plot_cos2 Logical. If `TRUE`, a plot of the cos2 values (quality of representation) is generated. Default is `FALSE`.
#' @param plot_corrplot Character. If not `NULL`, a correlation plot of the PCA results is generated. Options are `"cos2"`, `"contrib"`, or `"coord"`. Default is `NULL`.
#' @param plot_2d Logical. If `TRUE`, a 2D PCA plot is generated. Default is `TRUE`.
#' @param plot_3d Logical. If `TRUE`, a 3D PCA plot is generated. Default is `FALSE`.
#' @param gradient_cols A vector of colors for the gradient in the cos2 plot. Default is `c("#00AFBB", "#E7B800", "#FC4E07")`.
#'
#' @return A list containing:
#' \itemize{
#'   \item `pca_result`: The result of the `prcomp` function.
#'   \item `eigenvalues`: Eigenvalues of the principal components.
#'   \item `variance`: Variance explained by each principal component.
#'   \item `cumulative_variance`: Cumulative variance explained by the principal components.
#'   \item `cos2`: Cos2 values (quality of representation) for each sample.
#'   \item `ind`: The result of `get_pca_ind`, containing cos2, coord, and contrib.
#' }
#'
#' @examples
#' # Example 1: PCA analysis with a design matrix
#' # Example count matrix (samples as rows, genes as columns)
#' count_matrix <- data.frame(
#'   Gene1 = c(10, 20, 30),
#'   Gene2 = c(15, 25, 35),
#'   Gene3 = c(12, 22, 32)
#' )
#'
#' # Example design matrix (conditions for each sample)
#' design_matrix <- data.frame(
#'   Condition = c("A", "A", "B")
#' )
#'
#' # Perform PCA analysis and generate plots
#' pca_results <- pca_analysis(count_matrix, design_matrix = design_matrix, color_by = "Condition", shape_by = NULL, center = TRUE, scale = TRUE, plot_3d = TRUE)
#' 
#' # Access PCA results
#' pca_results$pca_result
#' pca_results$eigenvalues
#' pca_results$cos2
#'
#' # Example 2: PCA analysis with color and shape vectors (no design matrix)
#' color <- c("F", "F", "S")
#' shape <- c("FA", "FB", "SB")
#'
#' # Perform PCA analysis with custom color and shape vectors, generating 2D and 3D plots
#' pca_results <- pca_analysis(count_matrix, plot_2d = TRUE, plot_cos2 = FALSE, design_matrix = NULL, color_by = color, shape_by = shape, plot_3d = TRUE)

#' @importFrom factoextra fviz_pca_ind get_pca_ind get_eigenvalue
#' @importFrom plotly plot_ly
#' @importFrom corrplot corrplot
#' @importFrom stats prcomp
#' @importFrom grDevices colorRampPalette
#' @importFrom ggplot2 ggplot aes geom_point scale_color_manual scale_shape_manual theme_bw labs
#' @export
pca_analysis <- function(count_matrix, design_matrix = NULL, color_by = NULL, shape_by = NULL, center = TRUE, scale = FALSE, plot_cos2 = FALSE, plot_corrplot = NULL, plot_2d = TRUE, plot_3d = FALSE, gradient_cols = c("#00AFBB", "#E7B800", "#FC4E07")) {
  # Load necessary libraries
  require(factoextra)
  require(plotly)
  require(corrplot)
  require(ggplot2)
  
  # Ensure count_matrix is a data frame with samples as rows and genes as columns
  if (!is.data.frame(count_matrix)) {
    stop("Error: count_matrix must be a data frame with samples as row names and genes as columns.")
  }
  
  # Ensure row names are preserved
  if (is.null(rownames(count_matrix))) {
    stop("Error: count_matrix must have row names corresponding to sample names.")
  }
  
  # Check if design_matrix is provided and matches the number of samples
  if (!is.null(design_matrix)) {
    if (nrow(design_matrix) != nrow(count_matrix)) {
      stop("Error: design_matrix must have the same number of rows as count_matrix.")
    }
    # Ensure design_matrix is a data frame
    if (!is.data.frame(design_matrix)) {
      design_matrix <- as.data.frame(design_matrix)
    }
  }
  
  # Handle color_by
  if (!is.null(color_by)) {
    if (is.character(color_by) && length(color_by) == 1) {
      if (!is.null(design_matrix) && color_by %in% colnames(design_matrix)) {
        color_data <- design_matrix[[color_by]]
        color_title <- color_by
      } else {
        stop(paste("Error: color_by column", color_by, "not found in design_matrix."))
      }
    } else if (is.vector(color_by) && length(color_by) == nrow(count_matrix)) {
      color_data <- color_by
      color_title <- deparse(substitute(color_by))
    } else {
      stop("Error: color_by must be a column name in design_matrix or a vector of the same length as the number of samples.")
    }
  } else {
    color_data <- NULL
    color_title <- NULL
  }
  
  # Handle shape_by
  if (!is.null(shape_by)) {
    if (is.character(shape_by) && length(shape_by) == 1) {
      if (!is.null(design_matrix) && shape_by %in% colnames(design_matrix)) {
        shape_data <- design_matrix[[shape_by]]
        shape_title <- shape_by
      } else {
        stop(paste("Error: shape_by column", shape_by, "not found in design_matrix."))
      }
    } else if (is.vector(shape_by) && length(shape_by) == nrow(count_matrix)) {
      shape_data <- shape_by
      shape_title <- deparse(substitute(shape_by))
    } else {
      stop("Error: shape_by must be a column name in design_matrix or a vector of the same length as the number of samples.")
    }
  } else {
    shape_data <- NULL
    shape_title <- NULL
  }
  
  # Perform PCA
  pca_result <- prcomp(count_matrix, center = center, scale = scale)
  
  # Compute eigenvalues, variance, and cumulative variance
  eigenvalues <- get_eigenvalue(pca_result)
  variance <- eigenvalues$variance.percent
  cumulative_variance <- eigenvalues$cumulative.variance.percent
  
  # Get cos2 values (quality of representation)
  ind <- get_pca_ind(pca_result)
  cos2 <- ind$cos2
  
  # Plot cos2 values
  if (plot_cos2) {
    print(fviz_pca_ind(pca_result, 
                       col.ind = "cos2",  # Color by cos2 values
                       gradient.cols = gradient_cols,  # Gradient colors
                       repel = TRUE))  # Avoid text overlapping
  }
  
  # Plot correlation plot
  if (!is.null(plot_corrplot)) {
    if (plot_corrplot %in% c("cos2", "contrib", "coord")) {
      if (plot_corrplot == "cos2") {
        corrplot(ind$cos2, is.corr = FALSE, title = "Cos2 Values", tl.cex = 0.8, tl.col = "black", mar = c(0, 0, 1, 0))
      } else if (plot_corrplot == "contrib") {
        corrplot(ind$contrib, is.corr = FALSE, title = "Contributions", tl.cex = 0.8, tl.col = "black", mar = c(0, 0, 1, 0))
      } else if (plot_corrplot == "coord") {
        corrplot(ind$coord, is.corr = FALSE, title = "Coordinates", tl.cex = 0.8, tl.col = "black", mar = c(0, 0, 1, 0))
      }
    } else {
      stop("Invalid plot_corrplot option. Choose 'cos2', 'contrib', or 'coord'.")
    }
  }
  
  # Plot 2D PCA
  if (plot_2d) {
    pca_data <- as.data.frame(pca_result$x)
    pca_data$Sample <- rownames(pca_data)
    
    if (!is.null(color_data)) {
      pca_data$Color <- color_data
    }
    if (!is.null(shape_data)) {
      pca_data$Shape <- shape_data
    }
    
    # Create the 2D PCA plot
    pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2)) +
      geom_point(aes(color = if (!is.null(color_data)) Color else NULL,
                     shape = if (!is.null(shape_data)) Shape else NULL),
                 size = 3) +
      theme_bw() +
      labs(title = "PCA Plot", 
           x = paste0("PC1 (", round(variance[1], 2), "%)"), 
           y = paste0("PC2 (", round(variance[2], 2), "%)"))
    
    # Add custom colors if specified
    if (!is.null(color_data)) {
      n_colors <- length(unique(color_data))
      colors <- colorRampPalette(gradient_cols)(n_colors)
      pca_plot <- pca_plot + 
        scale_color_manual(values = colors, name = color_title)
    }
    
    # Add custom shapes if specified
    if (!is.null(shape_data)) {
      shapes <- rep(16:25, length.out = length(unique(shape_data)))
      pca_plot <- pca_plot + 
        scale_shape_manual(values = shapes, name = shape_title)
    }
    
    print(pca_plot)
  }
  
  # Plot 3D PCA
  if (plot_3d) {
    pca_3d_data <- as.data.frame(pca_result$x)
    
    if (!is.null(color_data)) {
      pca_3d_data$Color <- color_data
    }
    if (!is.null(shape_data)) {
      pca_3d_data$Shape <- shape_data
    }
    
    # Create a unique set of shape labels (e.g., FA, SA, UB) only if shape_by is provided
    if (!is.null(shape_data)) {
      unique_shapes <- unique(shape_data)
      
      # Map each unique shape to a plotly symbol
      valid_symbols <- c("circle", "square", "diamond", "cross", "x", "triangle-up", "triangle-down", "pentagon", "hexagon", "star")
      shape_symbol_map <- setNames(valid_symbols[1:length(unique_shapes)], unique_shapes)
      
      # Assign symbols to each shape in the data
      marker_symbols <- shape_symbol_map[as.character(shape_data)]
      
      # Define the color scale
      color_scale <- colorRampPalette(gradient_cols)(length(unique(color_data)))
      
      # Create the 3D plot
      pca_3d <- plot_ly(pca_3d_data, 
                        x = ~PC1, y = ~PC2, z = ~PC3,
                        color = ~Color,
                        colors = color_scale,
                        symbol = ~Shape,
                        symbols = valid_symbols[1:length(unique_shapes)],
                        text = rownames(pca_3d_data),  # Tooltip
                        type = 'scatter3d', 
                        mode = 'markers') 
      
      # Set the legend to use the custom shape labels
      pca_3d <- pca_3d %>% layout(
        legend = list(
          title = list(text = shape_title),
          itemsizing = 'constant',
          tracegroupgap = 10,
          traceorder = 'normal',
          font = list(size = 12)
        ),
        scene = list(
          xaxis = list(title = paste0("PC1 (", round(variance[1], 2), "%)")),
          yaxis = list(title = paste0("PC2 (", round(variance[2], 2), "%)")),
          zaxis = list(title = paste0("PC3 (", round(variance[3], 2), "%)"))
        )
      )
    } else {
      # If no shape data, just create the 3D plot without symbol
      pca_3d <- plot_ly(pca_3d_data, 
                        x = ~PC1, y = ~PC2, z = ~PC3,
                        color = ~Color,
                        colors = gradient_cols,
                        text = rownames(pca_3d_data),  # Tooltip
                        type = 'scatter3d', 
                        mode = 'markers') 
      
      pca_3d <- pca_3d %>% layout(
        scene = list(
          xaxis = list(title = paste0("PC1 (", round(variance[1], 2), "%)")),
          yaxis = list(title = paste0("PC2 (", round(variance[2], 2), "%)")),
          zaxis = list(title = paste0("PC3 (", round(variance[3], 2), "%)"))
        )
      )
    }
    
    # Display the plot
    print(pca_3d)
  }
  
  
  # Return results
  return(list(
    pca_result = pca_result,
    eigenvalues = eigenvalues,
    variance = variance,
    cumulative_variance = cumulative_variance,
    cos2 = cos2,
    ind = ind
  ))
}
