#' Compare Gene Association Rules Between Two Conditions
#'
#' This function compares two sets of gene association rules (as data frames),
#' identifies shared and unique rules, clusters them, and visualizes them as arc
#' or force-directed diagrams using `ggraph`. It also returns the rule data and 
#' cluster membership information for further analysis.
#'
#' @param rules_df1 A data frame of association rules from condition 1. Must contain columns `LHS` and `RHS`.
#' @param rules_df2 A data frame of association rules from condition 2. Must contain columns `LHS` and `RHS`.
#' @param condition_names A character vector of length 2 indicating the names of the two conditions. Defaults to `c("Condition1", "Condition2")`.
#' @param layout_type A character string indicating the layout for plotting. Use `"linear"` for arc diagrams (default) or `"fr"` for Fruchterman-Reingold force-directed layout.
#' @param comparison Compare complete gene sets (genes, default, as in the original relationship-focused approach), or complete expression-state itemsets (items). Metrics do not define rule identity. Original rule tables and gene-projection graphs are also returned.
#' @param plot Logical; if `TRUE`, plots the shared and unique rule graphs. Default is `TRUE`.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{shared_rules}{A data frame of rules common to both conditions.}
#'   \item{unique_rules_condition1}{Rules unique to condition 1.}
#'   \item{unique_rules_condition2}{Rules unique to condition 2.}
#'   \item{clusters}{A list of data frames containing gene-to-cluster mappings for shared and unique rules.}
#'   \item{plots}{A list of `ggplot` objects for shared and unique rule graphs.}
#' }
#'
#' @examples
#' \dontrun{
#' data <- matrix(rnorm(80), nrow=8, ncol=10)
#' colnames(data) <- paste0("Gene_", 1:10)
#' result <- mine_gene_associations(data, 1:4)
#' result_x <- mine_gene_associations(data, 5:8)
#' compare_gene_rules(rules_df1 = result$rules_df, rules_df2 = result_x$rules_df,
#'                    condition_names = c("F", "S"), plot = TRUE)
#' }
#'
#' @export
compare_gene_rules <- function(rules_df1, rules_df2, condition_names = c("Condition1", "Condition2"),
                               plot = TRUE, layout_type = "linear",
                               comparison = c("genes", "items")) {
  comparison <- match.arg(comparison)
  if (length(condition_names) != 2L || anyNA(condition_names)) stop("Supply two condition_names.")
  if (!is.logical(plot) || length(plot) != 1L || is.na(plot)) stop("plot must be TRUE or FALSE.")
  parse_side <- function(x) {
    x <- trimws(x)
    x <- sub("^\\{", "", x); x <- sub("\\}$", "", x)
    # Interval labels contain commas, so split only before a new gene=item token.
    if (!nzchar(x)) return(character())
    items <- trimws(strsplit(x, ",(?=\\s*[^,={}]+\\s*=)", perl = TRUE)[[1]])
    if (comparison == "genes") {
      if (!any(grepl("=", items, fixed = TRUE))) items <- trimws(strsplit(x, ",", fixed = TRUE)[[1]])
      items <- trimws(sub("=.*$", "", items))
    }
    sort(unique(items))
  }
  canonical <- function(df) {
    if (!is.data.frame(df)) stop("Rule inputs must be data frames.")
    if (!nrow(df)) return(data.frame(LHS = character(), RHS = character()))
    if (!all(c("LHS", "RHS") %in% names(df)) || anyNA(df$LHS) || anyNA(df$RHS)) stop("Rule tables require LHS and RHS without NA.")
    out <- data.frame(LHS = vapply(as.character(df$LHS), function(x) paste0("{", paste(parse_side(x), collapse = ","), "}"), character(1)),
      RHS = vapply(as.character(df$RHS), function(x) paste0("{", paste(parse_side(x), collapse = ","), "}"), character(1)))
    unique(out)
  }
  first <- canonical(rules_df1); second <- canonical(rules_df2)
  # Only the canonical itemsets define identity, never support/confidence/lift.
  keys <- function(df) paste(df$LHS, df$RHS, sep = "\r")
  common <- keys(first) %in% keys(second)
  shared <- first[common, , drop = FALSE]
  unique1 <- first[!common, , drop = FALSE]
  unique2 <- second[!keys(second) %in% keys(first), , drop = FALSE]
  process <- function(df) {
    if (!nrow(df)) return(list(plot = NULL, clusters = data.frame(Gene = character(), Cluster = integer()), graph = NULL))
    # Project every antecedent/consequent member to the gene graph; keep the
    # complete original rule separately so this projection does not erase it.
    edges <- do.call(rbind, lapply(seq_len(nrow(df)), function(i) {
      left <- parse_side(df$LHS[i]); right <- parse_side(df$RHS[i])
      left <- unique(sub("=.*$", "", left)); right <- unique(sub("=.*$", "", right))
      if (!length(left) || !length(right)) return(NULL)
      expand.grid(from = left, to = right, stringsAsFactors = FALSE)
    }))
    if (is.null(edges) || !nrow(edges)) return(list(plot = NULL, clusters = data.frame(Gene = character(), Cluster = integer()), graph = NULL))
    graph <- igraph::graph_from_data_frame(unique(edges), directed = TRUE)
    mem <- igraph::membership(igraph::cluster_walktrap(graph))
    clusters <- data.frame(Gene = names(mem), Cluster = as.integer(mem))
    graph <- igraph::set_vertex_attr(graph, "Cluster", value = as.integer(mem))
    p <- NULL
    if (plot) {
      colors <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(max(mem))
      p <- ggraph::ggraph(graph, layout = layout_type)
      p <- p + if (layout_type == "linear") ggraph::geom_edge_arc(alpha = .3) else ggraph::geom_edge_link(alpha = .3)
      p <- p + ggraph::geom_node_point(ggplot2::aes(color = factor(.data$Cluster)), size = 3) +
        ggraph::geom_node_text(ggplot2::aes(label = .data$name), repel = TRUE, size = 3) +
        ggplot2::scale_color_manual(values = colors) + ggplot2::theme_void() +
        ggplot2::theme(legend.position = "none")
    }
    list(plot = p, clusters = clusters, graph = graph)
  }
  results <- list(shared = process(shared), condition1 = process(unique1), condition2 = process(unique2))
  if (plot) {
    titles <- c(paste("Shared Rules -", paste(condition_names, collapse = " vs ")),
                paste("Unique Rules -", condition_names[1]), paste("Unique Rules -", condition_names[2]))
    for (i in seq_along(results)) if (!is.null(results[[i]]$plot)) print(results[[i]]$plot + ggplot2::ggtitle(titles[i]))
  }
  list(shared_rules = shared, unique_rules_condition1 = unique1, unique_rules_condition2 = unique2,
    clusters = lapply(results, `[[`, "clusters"),
    plots = list(shared_plot = results$shared$plot, unique_plot_condition1 = results$condition1$plot,
                 unique_plot_condition2 = results$condition2$plot),
    graphs = lapply(results, `[[`, "graph"), comparison = comparison,
    original_rules = list(condition1 = rules_df1, condition2 = rules_df2))
}
