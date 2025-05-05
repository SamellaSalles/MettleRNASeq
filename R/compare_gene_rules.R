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
#' @import igraph
#' @import dplyr
#' @import ggraph
#' @import RColorBrewer
#' @import ggplot2
#' @export
compare_gene_rules <- function(rules_df1, rules_df2,
                               condition_names = c("Condition1", "Condition2"),
                               plot = TRUE) {
  # Load required libraries
  library(igraph)
  library(dplyr)
  library(RColorBrewer)
  library(ggraph)
  library(ggplot2)
  
  # Clean LHS and RHS
  rules_df1 <- rules_df1 %>%
    mutate(LHS = gsub("\\=.*", "", LHS),
           RHS = gsub("\\=.*", "", RHS))
  rules_df2 <- rules_df2 %>%
    mutate(LHS = gsub("\\=.*", "", LHS),
           RHS = gsub("\\=.*", "", RHS))
  
  # Compare rules
  shared_rules <- intersect(rules_df1, rules_df2)[, c("LHS", "RHS")]
  unique_rules1 <- setdiff(rules_df1, rules_df2)[, c("LHS", "RHS")]
  unique_rules2 <- setdiff(rules_df2, rules_df1)[, c("LHS", "RHS")]
  
  # Helper to plot and extract cluster info
  process_rule_graph <- function(df) {
    if (nrow(df) == 0) return(list(plot = NULL, clusters = NULL))
    
    g <- graph_from_data_frame(df, directed = TRUE)
    com <- cluster_walktrap(g)
    mem <- membership(com)
    
    mycolor <- brewer.pal(n = max(mem), name = "Paired")
    V(g)$color <- mycolor[mem]
    
    plot <- ggraph(g, layout = "linear") +
      geom_edge_arc(edge_colour = "black", edge_alpha = 0.2, edge_width = 0.3, fold = TRUE) +
      geom_node_point(aes(size = 5, color = as.factor(color), fill = color), alpha = 0.5) +
      scale_size_continuous(range = c(0.5, 8)) +
      scale_color_manual(values = mycolor) +
      geom_node_text(aes(label = name), angle = 65, hjust = 1, nudge_y = -0.3, size = 2.3) +
      theme_void() +
      theme(legend.position = "none",
            plot.margin = unit(c(0, 0, 0.4, 0), "null"),
            panel.spacing = unit(c(0, 0, 3.4, 0), "null")) +
      expand_limits(x = c(-1.2, 1.2), y = c(-5.6, 1.2))
    
    list(plot = plot, clusters = data.frame(Gene = names(mem), Cluster = mem))
  }
  
  # Generate plots and clusters
  shared_res <- process_rule_graph(shared_rules)
  unique1_res <- process_rule_graph(unique_rules1)
  unique2_res <- process_rule_graph(unique_rules2)
  
  # Optionally plot
  if (plot) {
    print(shared_res$plot + ggtitle(paste("Shared Rules -", paste(condition_names, collapse = " vs "))))
    print(unique1_res$plot + ggtitle(paste("Unique Rules -", condition_names[1])))
    print(unique2_res$plot + ggtitle(paste("Unique Rules -", condition_names[2])))
  }
  
  # Return full result
  list(
    shared_rules = shared_rules,
    unique_rules_condition1 = unique_rules1,
    unique_rules_condition2 = unique_rules2,
    clusters = list(
      shared = shared_res$clusters,
      condition1 = unique1_res$clusters,
      condition2 = unique2_res$clusters
    ),
    plots = list(
      shared_plot = shared_res$plot,
      unique_plot_condition1 = unique1_res$plot,
      unique_plot_condition2 = unique2_res$plot
    )
  )
}
