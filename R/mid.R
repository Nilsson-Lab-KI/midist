#
# This file contains functions for MID correction, convolution,
# and related operations on MIDs
#


#' Isotopic enrichment of an MID (fraction of heavy atoms)
#'
#' @param mid and MID vector
#' @returns the enrichment value
isotopic_enrichment <- function(mid)
{
  n <- length(mid) - 1
  return(sum(mid * c(0:n)) / n)
}


#
# return a zero vector if the mid enrichment is below a threshold
#
filter_enrichment <- function(mids, tol = 0.01)
{
  if (isotopic_enrichment(mids) <= tol) {
    return(rep(0, length(mids)))
  }
  else {
    return(mids)
  }
}


#' Correct an MID vector for naturally occurring isotopes
# 13C correct
c13correct <- function(mid, constraint = TRUE)
{
  p <- 0.0107
  nrCarbon <- length(mid)-1
  end <- length(mid)
  correct <- matrix(0, end, end)

  for (d2 in 1:end){
    b1 <- c(0:(end-d2))
    b2 <- end-d2
    correct[d2:end, d2] <- binomvals(b1, b2, p)
  }

  return(pnnls(a = correct, b = mid)$x)
}

#' Geneerate a convolution matrix for an MID
#'
#' Convolution of two MIDs x * y can be written as a matrix multiplicationA(x).y
#' This function creates thematrix A(x) for a given number of carbon in y
#' @param x an MID
#' @param y_carbons the number of carbons in the vector y

# this creates the convolution matrix of size
# length(longer.mid) x length(middle.length.mid)
#convolution_matrix <- function(longer.mid, shorter.mid)

convolution_matrix <- function(x, y_carbons)
{
  #middle.length <- (length(longer.mid)-1) - (length(shorter.mid)-1) + 1
  x_carbons <- length(x) - 1
  # create empty matrix
  A <- matrix(0, x_carbons + y_carbons + 1, y_carbons + 1)
  # fill in columns
  for (i in 1:(y_carbons + 1)) {
    A[i:(i+x_carbons), i] <- x
  }
  return(A)
}

# find the non-negative least-squares solution to A y = z
# such that sum(y) = 1, y >= 0
nnls_solution <- function(A, z)
{
  return(pnnls(a = A, b = z, k = 0, sum = 1)$x)
}

# returns the solution y for x*y = z,
# where x and z are the shorter and longer mids, respectively
solution <- function(longer.mid, shorter.mid)
{
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

find_convolution <- function(longer.mid, shorter.mid, tol = 0.0107)
{
  A <- convolution_matrix(longer.mid, shorter.mid)
  sol <- nnls_solution(A, longer.mid)

  if(isotopic_enrichment(sol) < tol)
    return(NA)
  con <- A %*% sol
  if(isotopic_enrichment(con) < tol)
    return(NA)

  return(con)
}


