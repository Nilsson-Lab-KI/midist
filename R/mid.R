#
# This file contains functions for MID correction, convolution,
# and related operations on MIDs
#


#' Isotopic enrichment of an MID (fraction of heavy atoms)
#'
#' @param mid and MID vector
#' @returns the enrichment value
#' @export
isotopic_enrichment <- function(mid) {
  # number of carbon atoms
  n <- length(mid) - 1
  # compute and return isotopic enrichment
  return(sum(mid * c(0:n)) / n)
}


#
# return a zero vector if the mid enrichment is below a threshold
#
filter_enrichment <- function(mid, tol = 0.0107) {
  if (isotopic_enrichment(mid) <= tol) {
    return(rep(0, length(mid)))
  } else {
    return(mid)
  }
}

# calculate binom values
binomvals <- function(is, n, p) {
  bins <- c()
  for (j in 1:length(is)) {
    bin <- dbinom(is[j], n, p)
    bins <- c(bins, bin)
  }
  return(bins)
}

#' Correct an MID vector for naturally occurring isotopes
#' @param mid An MID vector to correct
#' @param p heavy atom natural abundance
#' @param constraint whether to constrain sum of corrected vector to 1
#' @export
c13correct <- function(mid, p = 0.0107, constraint = TRUE) {
  # number of carbon atoms
  nrCarbon <- length(mid) - 1
  # dimensions of the correction matrix
  end <- length(mid)
  # an empty correction matrix to be filled in by binom values
  correct <- matrix(0, end, end)

  # column-wise filling in the correction matrix
  for (d2 in 1:end) {
    b1 <- c(0:(end - d2))
    b2 <- end - d2
    correct[d2:end, d2] <- binomvals(b1, b2, p)
  }

  # if we do not have a constraint on the sum of the values
  if (constraint == FALSE) {
    return(pnnls(a = correct, b = mid)$x)
  }

  # if we have a sum 1 constraint (default)
  else {
    return(pnnls(a = correct, b = mid, sum = 1)$x)
  }
}


#' Correct an MID vector for naturally occurring isotopes and remove M+0
#' @param mid An MID vector to correct
#' @export
correct_and_remove <- function(mid) {
  return(c13correct(mid)[-1])
}


#' Add multiplicative gamma noise to a given MID from a given standard deviation
#' @param mid An MID vector
#' @param stdev standard deviation for the gamma distribution
#' Returns a noisy MID
#' @export
add_gamma_noise <- function(mid, stdev){
  # mean of the MID
  m <- mean(mid)
  # determine theta from stdev^2 = m * theta --> theta = stdev^2 / m
  theta <- (stdev^2) / m
  # determine k from stdev^2 = k * theta --> k = stdev^2 / theta
  k <- (stdev^2) / theta
  # generate a gamma distribution from these parameters
  gamma_vals <- rgamma(length(mid), k, scale = theta)
  # multiply the gamma values by the MID to make the noise multiplicative
  gamma_noise <- mid*gamma_vals
  noisy_mid <- (mid + gamma_noise) / sum(mid + gamma_noise)
  # stopifnot(sum(noisy_mid) == 1)
  return(noisy_mid)
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
convolute <- function(x, y) {
  return(as.vector(convolution_matrix(x, length(y) - 1) %*% y))
}


#' Convolute x*y for all column vectors y in matrix y_mat
#' @param x an MID vector
#' @param y_mat a matrix whose columns are MID vectors
#' @returns the matrix of convolution vectors x*y for each column y in y_mat
#' @export
convolute_cols <- function(x, y_mat) {
  return(convolution_matrix(x, nrow(y_mat) - 1) %*% y_mat)
}


# find the non-negative least-squares solution to A y = z
# such that sum(y) = 1, y >= 0
nnls_solution <- function(A, z) {
  return(pnnls(a = A, b = z, k = 0, sum = 1)$x)
}

# returns the solution y for x*y = z,
# where x and z are the shorter and longer mids, respectively
solution <- function(longer.mid, shorter.mid) {
  # RN: here there should probably be a check that
  # length(longer.mid) > length(shorter.mid)

  # DS: it will produce a matrix dimension error if length(longer.mid) < length(shorter.mid)
  A <- convolution_matrix(longer.mid, shorter.mid)

  # calculating the solution that minimizes the difference between the two MIDs
  # when the smaller metabolite convolutes to the larger one
  return(nnls_solution(A, longer.mid))
}

#' Find the optimal convolution for two MIDs
#'
#' Computes the convolution x*y (for x*y = z), where x and z are the shorter and
#' longer mids, respectively and y is the unknown mid
#' @param longer.mid the longer mid z
#' @param shorter.mid the shorter mid x
#' @param tol a threshold below which the convolution is not computed
#' @returns the convoluted MID vector, or NA if isotopoic enrichment is less than than tolerance
#' @export

find_convolution <- function(longer.mid, shorter.mid, tol = 0.0107) {
  A <- convolution_matrix(longer.mid, shorter.mid)
  sol <- nnls_solution(A, longer.mid)

  if (isotopic_enrichment(sol) < tol) {
    return(NA)
  }
  con <- A %*% sol
  if (isotopic_enrichment(con) < tol) {
    return(NA)
  }

  return(con)
}
