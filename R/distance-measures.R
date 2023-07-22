#
# distance and similarity measures
#


#' Cosine similarity between two vectors of the same length.
#' @param x a vector
#' @param y a vector
#' @returns the cosine similarity, or NA if either x or y is a zero vector
#' @export
cosine_sim <- function(x, y) {
  # compute norms
  x_norm <- sqrt(sum(x^2))
  y_norm <- sqrt(sum(y^2))
  # test for zero vector
  if (x_norm == 0 || y_norm == 0) {
    return(as.double(NA))
  } else {
    return(sum(x * y) / (x_norm * y_norm))
  }
}


#' Cosine similarity between two vectors of the same length without M+0.
#' @param x a vector
#' @param y a vector
#' @returns the cosine similarity, or NA if either x or y is a zero vector
#' @export
#'
cosine_sim_no_m0 <- function(x, y) apply_no_m0(cosine_sim, x, y)

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
  stopifnot(is.vector(x))
  stopifnot(is.vector(y))
  stopifnot(length(x) == length(y))
  return(1 - cosine_sim(x,y))
}


#' Cosine distance between two vectors of the same length without M+0.
#'
#' This is 1 minus the cosine similarity.
#'
#' @param x a vector
#' @param y a vector
#' @returns the cosine similarity, or NA if either x or y is a zero vector
#' @export
#'
cosine_dist_no_m0 <- function(x, y) apply_no_m0(cosine_dist, x, y)


#' Dot product similarity
#' @param x a vector
#' @param y a vector
#' @returns the dot (scalar) product x.y
#' @export
#'
dot_sim <- function(x, y) {
  return(sum(x * y))
}


#' Dot product similarity without M+0
#' @param x a vector
#' @param y a vector
#' @returns the dot (scalar) product x.y
#' @export
#'
dot_sim_no_m0 <- function(x, y) apply_no_m0(dot_sim, x, y)


#' Dot product distance (not really a distance)
#' @param x a vector
#' @param y a vector
#' @returns a distance based on the dot (scalar) product x.y
#' This is always nonnegative, and zero if x = y
#' @export
#'
dot_dist <- function(x, y)
{
  stopifnot(is.vector(x))
  stopifnot(is.vector(y))
  stopifnot(length(x) == length(y))
  return(sqrt(sum(x*x)*sum(y*y)) - sum(x*y))
}


#' Dot product distance (not really a distance) without M+0
#' @param x a vector
#' @param y a vector
#' @returns a distance based on the dot (scalar) product x.y
#' This is always nonnegative, and zero if x = y
#' @export
#'
dot_dist_no_m0 <- function(x, y) apply_no_m0(dot_dist, x, y)


#' The squared Euclidean distance (sum of squares)
#'
#' @param x a vector
#' @param y a vector
#' @returns the square of the Euclidean distance between x and y
#' @export
#'
euclidean_dist_sq <- function(x, y)
{
  stopifnot(is.vector(x))
  stopifnot(is.vector(y))
  stopifnot(length(x) == length(y))
  diff <- x - y
  return(sum(diff * diff))
}


#' The squared Euclidean distance (sum of squares) without M+0
#'
#' @param x a vector
#' @param y a vector
#' @returns the square of the Euclidean distance between x and y
#' @export
#'
euclidean_dist_sq_no_m0 <- function(x, y) apply_no_m0(euclidean_dist_sq, x, y)


#' The Euclidean distance
#'
#' @param x a vector
#' @param y a vector
#' @returns the Euclidean distance between x and y
#' @export
#'
euclidean_dist <- function(x, y) {
  return(sqrt(euclidean_dist_sq(x, y)))
}


#' The Euclidean distance without M+0
#'
#' @param x a vector
#' @param y a vector
#' @returns the Euclidean distance between x and y without M+0
#' @export
#'
euclidean_dist_no_m0 <- function(x, y) apply_no_m0(euclidean_dist, x, y)


#' The Jensen-Shannon (JS) distance
#'
#' @param x a vector
#' @param y a vector
#' @returns the JS distance between x and y
#' @export
#'
jensen_shannon <- function(x, y) {
  # Compute the average of the input vectors
  m <- (x + y) / 2

  # Compute the Kullback-Leibler divergences
  kl_x <- sum(x * log2(x / m))
  kl_y <- sum(y * log2(y / m))

  # Compute the Jensen-Shannon distance
  jsd <- sqrt((kl_x + kl_y) / 2)

  return(jsd)
}

#' @export
manhattan_distance <- function(x, y) {
  return(sum(abs(x - y)))
}

#' @export
canberra_distance <- function(x, y) {
  return(sum(abs(x - y) / (abs(x) + abs(y))))
}

#' @export
bray_curtis_distance <- function(x, y) {
  return(sum(abs(x - y)) / sum(abs(x + y)))
}


#' @export
mi_weighted_distance <- function(x, y, p = 1) {
  n <- length(x) - 1
  return((sum(abs(x - y)^p * c(0:n)))^(1/p))
}


#' @export
mi_weighted_dist_normalized <- function(x, y, p = 1) {
  n <- length(x) - 1
  return(((1/n)*sum(abs(x - y)^p * c(0:n)))^(1/p))
}
