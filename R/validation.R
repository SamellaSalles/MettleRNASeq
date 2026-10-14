.mettle_integer <- function(x, name, minimum = 1L) {
  if (length(x) != 1L || !is.numeric(x) || is.na(x) || !is.finite(x) ||
      x != floor(x) || x < minimum) {
    stop(name, " must be a single integer >= ", minimum, ".", call. = FALSE)
  }
  as.integer(x)
}

.mettle_numeric_matrix <- function(x, name = "data", nonnegative = FALSE) {
  if (!(is.matrix(x) || is.data.frame(x)) || length(dim(x)) != 2L ||
      nrow(x) == 0L || ncol(x) == 0L) {
    stop(name, " must be a non-empty numeric matrix or data frame.", call. = FALSE)
  }
  if (is.data.frame(x) && !all(vapply(x, is.numeric, logical(1)))) {
    stop(name, " must contain numeric columns only.", call. = FALSE)
  }
  x <- as.matrix(x)
  if (!is.numeric(x) || any(!is.finite(x))) {
    stop(name, " must contain finite numeric values without NA.", call. = FALSE)
  }
  if (nonnegative && any(x < 0)) stop(name, " must be non-negative.", call. = FALSE)
  x
}

.mettle_names <- function(x, name) {
  if (is.null(x) || anyNA(x) || any(!nzchar(x)) || anyDuplicated(x)) {
    stop(name, " must be present, non-empty and unique.", call. = FALSE)
  }
  invisible(x)
}

.mettle_split <- function(data, condition_split) {
  if (!is.numeric(condition_split) || length(condition_split) != 2L ||
      anyNA(condition_split) || any(!is.finite(condition_split)) ||
      any(condition_split != floor(condition_split)) ||
      condition_split[1] < 2L || condition_split[2] != nrow(data) ||
      condition_split[2] - condition_split[1] < 2L) {
    stop("condition_split must define two groups of at least two samples and end at nrow(data).",
         call. = FALSE)
  }
  list(seq_len(condition_split[1]), seq.int(condition_split[1] + 1L, condition_split[2]))
}
