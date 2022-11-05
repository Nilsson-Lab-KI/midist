#
# distance and similarity measures
#

#' Cosine similarity between two vectors of the same length.
#' @param x a vector
#' @param y a vector
#' @returns the cosine similarity, or NA if either x or y is a zero vector
#' @export
cosine_sim <- function(x, y)
{
  cosi <- sum(x * y) / (sqrt(sum(x^2)) * sqrt(sum(y^2)))
  return(cosi)
}

#' Cosine distance between two vectors of the same length.
#'
#' This is 1 minus the cosine similarity.
#'
#' @param x a vector
#' @param y a vector
#' @returns the cosine similarity, or NA if either x or y is a zero vector
#' @export
cosine_dist <- function(x, y)
{
  return(1 - cosine_sim)
}

#' Apply a function on MIDs x,y after removing M+0
apply_no_m0 <- function(f, x, y)
{
  return(f(x[-1], y[-1]))
}

calc_dot <- function(x, y) {1 - sum(x*y)}
