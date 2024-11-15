% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convolution-distance.R
\name{filter_pairwise_matrix_global}
\alias{filter_pairwise_matrix_global}
\title{Filter a matrix so that a given fraction of its real-valued elements
are set to the corresponding quantile}
\usage{
filter_pairwise_matrix_global(pairwise_matrix, percentile = 0.01)
}
\arguments{
\item{pairwise_matrix}{A symmetric matrix}

\item{percentile}{The fraction (NOT a percentile) of elements to keep}
}
\description{
Filter a matrix so that a given fraction of its real-valued elements
are set to the corresponding quantile
}
