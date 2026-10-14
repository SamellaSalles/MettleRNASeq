#' Inspect Gene Coverage and Recurrence Across Resamples
#'
#' Summarize which genes were drawn, successfully evaluated and present in resamples
#' meeting a score threshold. Counts are per unique gene and resample, so repeated
#' draws of a gene within one resample do not inflate its recurrence.
#'
#' @param classification_result The complete result from classification_approach_advanced().
#' @param cutoff Optional minimum score. With NULL, qualification counts and fractions are NA.
#' @param metric_column Resample summary column, default summary_test_accuracy.
#' @param algorithm Optional algorithm names; NULL summarizes all algorithms separately.
#' @param gene_universe Optional complete gene-ID vector. By default, use the gene_order
#' metadata retained on genes_table. Supply this vector if that metadata was lost during export.
#' @return A list containing gene_table, coverage_table, resample_table and settings.
#' gene_table includes all genes, even those never sampled. qualification_fraction uses
#' the number of successfully evaluated resamples containing each gene as its denominator.
#' qualified_presence_fraction uses the number of qualifying resamples. An undefined
#' denominator produces NA. The tables report descriptive recurrence, not gene-level
#' importance, significance, independent validation or calibrated selection stability.
#' @examples
#' \dontrun{
#' diagnostics <- summarize_gene_resampling(classification, cutoff = 0.85)
#' diagnostics$coverage_table
#' diagnostics$gene_table
#' }
#' @export
summarize_gene_resampling <- function(classification_result, cutoff = NULL,
                                     metric_column = "summary_test_accuracy", algorithm = NULL,
                                     gene_universe = NULL) {
  if (!is.list(classification_result) ||
      !is.data.frame(classification_result$eval_table) ||
      !is.data.frame(classification_result$gene_ids_table)) {
    stop("Supply a complete classification result with eval_table and gene_ids_table.")
  }
  if (!is.character(metric_column) || length(metric_column) != 1L ||
      is.na(metric_column) || !nzchar(metric_column)) stop("metric_column must name one resample summary column.")
  metric_column <- tolower(metric_column)
  values <- classification_result$eval_table
  names(values) <- tolower(names(values))
  if (!nrow(values) || !all(c("algorithm", "resample", metric_column) %in% names(values))) {
    stop("Evaluation data must contain algorithm, resample and the requested summary column.")
  }
  if (!is.numeric(values[[metric_column]]) || anyNA(values$algorithm) ||
      any(!nzchar(as.character(values$algorithm)))) stop("Invalid algorithm labels or summary values.")
  if (!is.null(cutoff) && (!is.numeric(cutoff) || length(cutoff) != 1L ||
      !is.finite(cutoff))) stop("cutoff must be NULL or one finite numeric value.")
  methods <- unique(as.character(values$algorithm))
  if (!is.null(algorithm)) {
    if (!is.character(algorithm) || !length(algorithm) || anyNA(algorithm) ||
        anyDuplicated(algorithm) || !all(algorithm %in% methods)) stop("Unknown or duplicated algorithm.")
    methods <- algorithm
  }
  if (is.null(gene_universe)) gene_universe <- attr(classification_result$genes_table, "gene_order")
  if (is.null(gene_universe) || !length(gene_universe)) {
    stop("Supply the complete gene_universe when gene_order metadata is unavailable.")
  }
  gene_universe <- as.character(gene_universe)
  .mettle_names(gene_universe, "Gene universe")
  draws <- classification_result$gene_ids_table
  .mettle_names(names(draws), "Resample names")
  if (!nrow(draws) || !ncol(draws)) stop("Gene draws must be non-empty.")
  gene_sets <- lapply(draws, function(ids) {
    ids <- as.character(ids)
    if (anyNA(ids) || any(!ids %in% gene_universe)) stop("Gene draws must belong to the complete gene universe.")
    unique(ids)
  })
  gene_rows <- list(); coverage_rows <- list(); resample_rows <- list()
  for (method in methods) {
    table <- values[as.character(values$algorithm) == method, , drop = FALSE]
    ids <- unique(table$resample)
    if (!is.numeric(ids) || anyNA(ids) || any(!is.finite(ids)) ||
        any(ids < 1 | ids != floor(ids))) stop("Resample IDs must be positive integers.")
    ids <- sort(ids)
    columns <- paste0("Resample_", ids)
    if (!all(columns %in% names(gene_sets))) stop("Evaluation resamples do not match gene draws.")
    scores <- vapply(ids, function(id) {
      rows <- table[table$resample == id, , drop = FALSE]
      unique_scores <- unique(rows[[metric_column]])
      if (length(unique_scores) != 1L) stop("The requested metric must have one consistent summary per resample.")
      if ("status" %in% names(rows) &&
          (anyNA(rows$status) || any(rows$status != "ok"))) return(NA_real_)
      if (!is.finite(unique_scores)) return(NA_real_)
      unique_scores
    }, numeric(1))
    evaluated <- is.finite(scores)
    qualified <- if (is.null(cutoff)) rep(NA, length(ids)) else evaluated & scores >= cutoff
    count_genes <- function(keep) {
      as.integer(tabulate(match(unlist(gene_sets[columns[keep]], use.names = FALSE),
        gene_universe), nbins = length(gene_universe)))
    }
    sampled_counts <- count_genes(rep(TRUE, length(ids)))
    evaluated_counts <- count_genes(evaluated)
    qualified_counts <- if (is.null(cutoff)) rep(NA_integer_, length(gene_universe)) else count_genes(qualified)
    fractions <- rep(NA_real_, length(gene_universe))
    if (!is.null(cutoff)) fractions[evaluated_counts > 0L] <-
      qualified_counts[evaluated_counts > 0L] / evaluated_counts[evaluated_counts > 0L]
    n_qualified <- if (is.null(cutoff)) NA_integer_ else sum(qualified)
    presence <- if (!is.na(n_qualified) && n_qualified > 0L)
      qualified_counts / n_qualified else rep(NA_real_, length(gene_universe))
    gene_rows[[method]] <- data.frame(algorithm = method, gene_id = gene_universe,
      resamples_sampled = sampled_counts, resamples_evaluated = evaluated_counts,
      resamples_qualified = qualified_counts, sampling_fraction = sampled_counts / length(ids),
      qualification_fraction = fractions, qualified_presence_fraction = presence)
    coverage_rows[[method]] <- data.frame(algorithm = method,
      total_resamples = length(ids), evaluated_resamples = sum(evaluated),
      qualifying_resamples = n_qualified, gene_universe_size = length(gene_universe),
      sampled_genes = sum(sampled_counts > 0L), sampling_coverage = mean(sampled_counts > 0L),
      evaluated_genes = sum(evaluated_counts > 0L), evaluated_coverage = mean(evaluated_counts > 0L),
      min_resamples_per_gene = min(sampled_counts), max_resamples_per_gene = max(sampled_counts))
    resample_rows[[method]] <- data.frame(algorithm = method, resample = ids,
      score = scores, evaluated = evaluated, qualified = qualified,
      unique_genes = lengths(gene_sets[columns]))
  }
  list(gene_table = do.call(rbind, gene_rows), coverage_table = do.call(rbind, coverage_rows),
    resample_table = do.call(rbind, resample_rows), settings = list(cutoff = cutoff,
      metric_column = metric_column, algorithms = methods, gene_universe_size = length(gene_universe)))
}
