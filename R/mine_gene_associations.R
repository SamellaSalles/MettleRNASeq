#' Mine Gene Associations with Rule-Based Algorithms and Visualizations
#'
#' This function performs association rule mining on gene expression data using either the Apriori or FP-Growth algorithm. 
#' It provides filtering options (support, confidence, lift, rule length) and multiple visualization styles.
#'
#' @param data A data frame or matrix of gene expression (samples in rows, genes in columns).
#' @param condition_rows A numeric vector specifying which rows (samples) to include in the analysis.
#' @param method The mining algorithm to use. Options are `"apriori"` (default) or `"fpgrowth"`.
#' @param supp Minimum support threshold (default: 0.5).
#' @param conf Minimum confidence threshold (default: 0.8).
#' @param lift Optional minimum lift threshold to filter rules (default: `NULL`, meaning no lift filtering).
#' @param minlen Minimum length of rules (default: 2).
#' @param maxlen Maximum length of rules (default: 10).
#' @param target_genes Optional character vector of gene names to focus on.
#' @param search_type Where to match `target_genes`: `"both"` (default), `"lhs"` (antecedent), or `"rhs"` (consequent).
#' @param plot_type The type of plot: `"scatter"`, `"graph"`, `"paracoord"`, or `"grouped"` (default: `"scatter"`).
#' @param interactive Logical. Use interactive visualizations if `TRUE` (default).
#' @param top_n The number of top rules to include in the plot (default: 100).
#' @param measure Plotting measure for `grouped` plots (`"support"`, `"confidence"`, `"lift"`, or `"order"`). Default is `"support"`.
#' @param shading Shading measure for `grouped` plots (default: `"lift"`).
#' @param color Base color for `grouped` plots (default: `"#0073e6"`).
#' @param lhs_max Maximum items to show on the left-hand side (LHS) in grouped plots (default: 50).
#' @param rhs_max Maximum items to show on the right-hand side (RHS) in grouped plots (default: 50).
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{transactions}{The `transactions` object used for rule mining.}
#'   \item{rules}{The mined association rules (of class `rules`).}
#'   \item{rules_df}{A data frame representation of the rules.}
#'   \item{plot}{The generated plot object, if applicable.}
#' }
#'
#' @details
#' The function automatically installs required packages if they are not already installed. The `lift` filter 
#' is applied after computing the `lift` measure explicitly to ensure compatibility across methods.
#'
#' For `"fpgrowth"`, rules are mined using the `fim4r()` function from the `arules` package.
#'
#' @examples
#' \dontrun{
#' data <- matrix(rnorm(80), nrow=8, ncol=10)
#' colnames(data) <- paste0("Gene_", 1:10)
#' result <- mine_gene_associations(data, 1:4)
#' result <- mine_gene_associations(data, 1:4, plot_type = "paracoord", top_n = 5)
#' result <- mine_gene_associations(data, 1:4, method = "fpgrowth", plot_type = "graph", interactive = FALSE, lift = 2)
#' 
#' #Focusing on genes of interest
#' result <- mine_gene_associations(data, 1:4, plot_type = "graph", interactive = TRUE, top_n = 10, target_genes = "Gene_1", search_type = "lhs")
#'}
#'
#' @import arules
#' @import arulesViz
#' @importFrom utils install.packages
#' @export
mine_gene_associations <- function(data, 
                                   condition_rows,
                                   method = c("apriori", "fpgrowth"),
                                   supp = 0.5,
                                   conf = 0.8,
                                   lift = NULL,
                                   minlen = 2,
                                   maxlen = 10,
                                   target_genes = NULL,
                                   search_type = c("both", "lhs", "rhs"),
                                   plot_type = c("scatter", "graph", "paracoord", "grouped"),
                                   interactive = TRUE,
                                   top_n = 100,
                                   measure = "support",
                                   shading = "lift",
                                   color = "#0073e6",
                                   lhs_max = 50,
                                   rhs_max = 50) {
  
  # Ensure required packages are installed and loaded
  required_pkgs <- c("arules", "arulesViz")
  for (pkg in required_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  }
  suppressPackageStartupMessages({
    library(arules)
    library(arulesViz)
  })
  
  # Validate inputs
  method <- match.arg(method)
  plot_type <- match.arg(plot_type)
  search_type <- match.arg(search_type)
  
  if (interactive && plot_type == "scatter") {
    if (!requireNamespace("plotly", quietly = TRUE)) install.packages("plotly")
  }
  if (interactive && plot_type == "graph") {
    if (!requireNamespace("htmlwidgets", quietly = TRUE)) install.packages("htmlwidgets")
  }
  if (!interactive && plot_type == "graph") {
    if (!requireNamespace("igraph", quietly = TRUE)) install.packages("igraph")
  }
  
  
  if (length(condition_rows) == 0) stop("condition_rows must specify at least one row")
  
  # Prepare data
  data <- as.data.frame(data)
  if (is.null(colnames(data))) colnames(data) <- paste0("Gene_", 1:ncol(data))
  
  cond_data <- data[condition_rows, ]
  transactions <- as(cond_data, "transactions")
  
  # Mine rules
  cat("Mining association rules using", method, "...\n")
  if (method == "apriori") {
    rules <- apriori(transactions, 
                     parameter = list(supp = supp, conf = conf,
                                      minlen = minlen, maxlen = maxlen,
                                      target = "rules"))
  } else if (method == "fpgrowth") {
    rules <- fim4r(transactions, 
                   method = "fpgrowth", 
                   target = "rules", 
                   supp = supp, 
                   conf = conf, 
                   zmin = minlen, 
                   zmax = maxlen)
  }
  
  # Filter by lift
  if (!is.null(lift)) {
    rules <- rules[quality(rules)$lift >= lift]
    cat("Filtered rules using lift ≥", lift, "\n")
  }
  
  # Focus on specific genes if requested
  if (!is.null(target_genes)) {
    target_items <- paste(target_genes, collapse = "|")
    rules <- switch(search_type,
                    "both" = subset(rules, lhs %pin% target_items | rhs %pin% target_items),
                    "lhs"  = subset(rules, lhs %pin% target_items),
                    "rhs"  = subset(rules, rhs %pin% target_items))
  }
  
  # Remove redundant rules
  if (length(rules) > 0) rules <- rules[!is.redundant(rules)]
  
  # Generate plot
  plot_obj <- NULL
  if (length(rules) > 0) {
    plot_rules <- if (length(rules) > top_n) head(sort(rules, by = "lift"), top_n) else rules
    
    cat("Preparing", plot_type, "plot...\n")
    plot_obj <- tryCatch({
      switch(plot_type,
             
             "scatter" = {
               engine <- if (interactive) "plotly" else "ggplot2"
               plot(plot_rules, method = "scatter", engine = engine, max = top_n)
             },
             
             "graph" = {
               engine <- if (interactive) "htmlwidget" else "igraph"
               plot(plot_rules, method = "graph", engine = engine, max = top_n)
             },
             
             "paracoord" = {
               plot(plot_rules, method = "paracoord", max = top_n)
             },
             
             "grouped" = {
               plot(plot_rules, 
                    method = "grouped", 
                    measure = measure,
                    shading = shading,
                    control = list(k = lhs_max, rhs_max = rhs_max),
                    col = color)
             }
      )
    }, error = function(e) {
      message("Plotting failed: ", e$message)
      NULL
    })
  }
  
  if (!is.null(plot_obj)) print(plot_obj)
  
  # Return results
  list(
    transactions = transactions,
    rules = rules,
    rules_df = if (length(rules) > 0) DATAFRAME(rules) else data.frame(),
    plot = plot_obj
  )
}
