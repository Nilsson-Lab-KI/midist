#
# This file contains functions for MID correction, convolution,
# and related operations on MIDs
#


#' Isotopic enrichment of an MID (fraction of heavy atoms)
#'
#' @param mid An MID vector
#' @returns The enrichment value
#' @export
isotopic_enrichment <- function(mid) {
  # number of carbon atoms
  n <- length(mid) - 1
  # compute and return isotopic enrichment
  return(sum(mid * c(0:n)) / n)
}


#' Return a zero vector if the mid enrichment is below a threshold
#' @param mid An MID vector
#' @param tol A threshold value
#' @returns The given mid vector, or a vector of zeros of the same length,
#'if the isotopic enrichment is less than threshold
filter_enrichment <- function(mid, tol = 0.0107)
{
  if (isotopic_enrichment(mid) <= tol) {
    return(rep(0, length(mid)))
  } else {
    return(mid)
  }
}


#' Correct an MID vector for naturally occurring isotopes
#'
#' @param mid An MID vector to correct
#' @param p The heavy atom natural abundance
#' @param constraint whether to constrain sum of corrected vector to 1
#' @export
c13correct <- function(mid, p = 0.0107, constraint = TRUE)
{
  if (!all(is.na(mid))){
    # dimensions of the correction matrix
    end <- length(mid)
    # an empty correction matrix to be filled in by binom values
    correct <- matrix(0, end, end)
    
    # column-wise filling in the correction matrix
    for (j in 1:end) {
      correct[j:end, j] <- dbinom(c(0:(end-j)), end-j, p)
    }
    
    # if we do not have a constraint on the sum of the values
    if (constraint == FALSE) {
      return(pnnls(a = correct, b = mid)$x)
    }
    # if we have a sum 1 constraint (default)
    else {
      return(pnnls(a = correct, b = mid, sum = 1)$x)
    } 
  } else return(mid)
  
}


#' Add multiplicative gamma noise to a given MID from a given standard deviation
#'
#' @param mid An MID vector
#' @param stdev standard deviation for the gamma distribution
#' Computes the gamma distribution parameters theta and k from the
#' mean of the MID and input standard deviation.
#' Multiplies the gamma values by the MID itself (multiplicative)
#' Adds this noise to the MID and returns the noisy MID
#' @export
add_gamma_noise <- function(mid, stdev)
{
  # mean of the MID
  m <- mean(mid)
  # determine theta from stdev^2 = m * theta --> theta = stdev^2 / m
  theta <- (stdev^2) / m
  # determine k from stdev^2 = k * theta --> k = stdev^2 / theta
  k <- (stdev^2) / (theta^2)
  # generate a gamma distribution from these parameters
  gamma_vals <- stats::rgamma(length(mid), k, scale = theta)
  # multiply the gamma values by the MID to make the noise multiplicative
  gamma_noise <- mid*gamma_vals
  noisy_mid <- (mid + gamma_noise) / sum(mid + gamma_noise)
  # stopifnot(sum(noisy_mid) == 1)
  return(noisy_mid)
}


#' Sample from a Dirichlet distribution.
#'
#' Draws n samples from a Dirichlet distribution with the given mean vector
#' and precision. Each component i of the resulting random vector has
#' marginal variance x_i*(1-x_i) / (1 + precision) where x_i is the i'th component
#' of the mean vector.
#'
#' @param mean Mean vector of the Dirichlet distribution. Must sum to 1.
#' @param precision Precision of the Dirichlet distribution,
#' @param n number of samples to draw
#' @returns A matrix with n columns, each an MID of the same length as mean
#'
random_dirichlet <- function(mean, precision, n)
{
  # If Y_i is distributed Gamma(shape = x_i*p, 1), i = 1, ..., n
  # then Y / sum(Y) is distributed Dirichlet(x_i*p, .. x_i*p),
  # with means E Y_i = x_i*p / sum(x_i*p) = x_i
  # and component variances Var Y_i = x_i(1 - x_i) / (1 + p)
  gamma_obs <- sapply(
    mean*precision,
    function(shape) stats::rgamma(n, shape = shape, scale = 1))
  return(t(gamma_obs / rowSums(gamma_obs)))
}


#' Sample MIDs from a Dirichlet distributon.
#'
#' Samples n MIDs from a Dirichlet distribution with the given mean vector
#' and a precision such that the standard deviation of an MI fraction
#' with (marginal) expectation = 1/2 equals stddev.
#'
#' @param mean An MID vector
#' @param stdev A standard deviation
#' @param n number of samples to draw
#' @returns A matrix with n columns, each an MID of the same length as mean
#' @export
random_mid <- function(mean, stdev, n)
{
  # If a component i of the random vector has mean 1/2, then
  # its variance is stdev^2 = 1/4 / (1 + p) where p is the precision
  # ==> p = 1/(4*stdev^2) - 1
  return(
    random_dirichlet(
      mean = mean,
      precision = 0.25 / (stdev^2) - 1,
      n))
}


#' Generate a convolution matrix for an MID
#'
#' Convolution of two MIDs x * y can be written as a matrix multiplication
#' A(x).y
#' This function creates the matrix A(x) for a given number of carbons in y
#' @param x an MID
#' @param y_carbons the number of carbons in the vector y
convolution_matrix <- function(x, y_carbons) {
  x_carbons <- length(x) - 1
  # create empty matrix
  A <- matrix(0, x_carbons + y_carbons + 1, y_carbons + 1)
  # fill in columns
  for (i in 1:(y_carbons + 1)) {
    A[i:(i + x_carbons), i] <- x
  }
  return(A)
}

#' Compute the convolution x*y of two MIDs x,y
#' @param x an MID vector
#' @param y an MID vector
#' @returns the convolution MID vector x*y
#' @export
convolute <- function(x, y)
{
  return(as.vector(convolution_matrix(x, length(y) - 1) %*% y))
}


#' Convolute x*y for all column vectors y in matrix y_mat
#' @param x an MID vector
#' @param y_mat a matrix whose columns are MID vectors
#' @returns the matrix of convolution vectors x*y for each column y in y_mat
convolute_cols <- function(x, y_mat) {
  return(convolution_matrix(x, nrow(y_mat) - 1) %*% y_mat)
}

#' Convolute columns of x with second dimension of y
convolute_array <- function(x_mat, y_array)
{
  # this yields an MI x experiments x z_index array
  n_atoms_x <- nrow(x_mat) - 1
  n_atoms_y <- dim(y_array)[1] - 1
  n_exp <- dim(y_array)[2]
  n_y <- dim(y_array)[3]
  return(
    aperm(
      vapply(
        1:n_exp,
        function(i) convolute_cols(x_mat[, i], y_array[, i, , drop = FALSE]),
        FUN.VALUE = array(0, dim = c(n_atoms_x + n_atoms_y + 1, n_y))
      ),
      c(1, 3, 2)
    )
  )
}


# find the non-negative least-squares solution to A y = z
# such that sum(y) = 1, y >= 0
nnls_solution <- function(A, z)
{
  return(pnnls(a = A, b = z, sum = 1)$x)
}


# returns the MID y minimizing || x*y = z ||
# where x and z are MIDs and z is longer than x
solution <- function(z, x)
{
  # number of atoms in MIDs
  n_x <- length(x) - 1
  n_z <- length(z) - 1
  # z must be larger than x
  stopifnot(n_x < n_z)
  n_y <- n_z - n_x
  A <- convolution_matrix(x, n_y)
  # solve the problem min || A.y - z || s.t. y >= 0 and sum(y) = 1
  return(nnls_solution(A, z))
}


#' Find the optimal convolution for two MIDs
#'
#' Computes the convolution x*y (for x*y = z), where x and z are the shorter and
#' longer MIDs, respectively, and y is the unknown MID
#' @param z the longer MID
#' @param x the shorter MID
#' @param tol a threshold below which the convolution is not computed
#' @returns the convoluted MID vector, or NA if isotopic enrichment is
#' less than than tolerance
#' @export
find_convolution <- function(z, x, tol = 0.0107)
{
  # this is redundant with solution() but avoids re-computing the matrix A ...
  n_x <- length(x) - 1
  n_z <- length(z) - 1
  # z must be larger than x
  stopifnot(n_x < n_z)
  n_y <- n_z - n_x
  A <- convolution_matrix(x, n_y)
  y <- nnls_solution(A, z)

  if(isotopic_enrichment(y) < tol)
    return(NA)
  # the convolution x*y
  con <- as.vector(A %*% y)
  if(isotopic_enrichment(con) < tol)
    return(NA)
  else
    return(con)
}
