% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prepare_input_matrix.R
\name{prepare_input_matrix}
\alias{prepare_input_matrix}
\title{Prepare Input Matrix for Classification}
\usage{
prepare_input_matrix(
  count_matrix,
  design_table = NULL,
  sample_col = NULL,
  condition_col = NULL,
  sample_vector = NULL,
  condition_vector = NULL
)
}
\arguments{
\item{count_matrix}{A matrix or data frame of count data, where rows are genes and columns are samples.}

\item{design_table}{A data frame containing sample and condition information. Default is NULL.}

\item{sample_col}{The name of the column in \code{design_table} containing sample names. Default is NULL.}

\item{condition_col}{The name of the column in \code{design_table} containing condition names. Default is NULL.}

\item{sample_vector}{A vector of sample names. Default is NULL.}

\item{condition_vector}{A vector of condition names corresponding to \code{sample_vector}. Default is NULL.}
}
\value{
A data frame in wide format with genes as rows, replicates as columns (named "Rep1", "Rep2", etc.),
and conditions as the last column.
}
\description{
This function prepares an input matrix for classification by organizing count data
into a wide format with consistent replicate naming.
}
\details{
This function:
\itemize{
\item Ensures all conditions have the same number of replicates
\item Converts the input matrix to a properly formatted data frame
\item Uses underscore-separated replicate names ("Rep_1" instead of "Replicate 1")
\item Maintains original gene order while arranging by condition
}
}
\examples{
# Example with vectors
sample <- c("F1", "F2", "F3", "F4", "S1", "S2", "S3", "S4")
condition <- c("F", "F", "F", "F", "S", "S", "S", "S")
result <- matrix(runif(8 * 100), nrow = 100, ncol = 8)
colnames(result) <- sample
prepare_input_matrix(count_matrix = result, sample_vector = sample, condition_vector = condition)

}
