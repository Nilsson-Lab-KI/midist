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
  # number of carbon atoms
  n <- length(mid) - 1
  # compute and return isotopic enrichment
  return(sum(mid * c(0:n)) / n)
}


#
# return a zero vector if the mid enrichment is below a threshold
#
filter_enrichment <- function(mid, tol = 0.0107)
{
  if (isotopic_enrichment(mid) <= tol) {
    return(rep(0, length(mid)))
  }
  else {
    return(mid)
  }
}

# calculate binom values
binomvals <- function(is, n, p){
  bins <- c()
  for (j in 1:length(is)){
    bin <- dbinom(is[j], n, p)
    bins <- c(bins, bin)
  }
  return(bins)
}

#' Correct an MID vector for naturally occurring isotopes
#' @param mid An MID vector to correct
#' @param constraint NOT CURRENTLY USED
#' @export
c13correct <- function(mid, p = 0.0107, constraint = TRUE)
{
  # number of carbon atoms
  nrCarbon <- length(mid)-1
  # dimensions of the correction matrix
  end <- length(mid)
  # an empty correction matrix to be filled in by binom values
  correct <- matrix(0, end, end)

  # column-wise filling in the correction matrix
  for (d2 in 1:end){
    b1 <- c(0:(end-d2))
    b2 <- end-d2
    correct[d2:end, d2] <- binomvals(b1, b2, p)
  }
  
  # if we do not have a constraint on the sum of the values
  if (constraint == FALSE){
    return(pnnls(a = correct, b = mid)$x)
  }
  
  # if we have a sum 1 constraint (default)
  else 
    return(pnnls(a = correct, b = mid, sum = 1)$x)

}


#' Correct an MID vector for naturally occurring isotopes and remove M+0
#' @param mid An MID vector to correct
#' @export
correct_and_remove <- function(mid){
  return(c13correct(mid)[-1])
}

#' Generate a convolution matrix for an MID
#'
#' Convolution of two MIDs x * y can be written as a matrix multiplication
#' A(x).y
#' This function creates the matrix A(x) for a given number of carbons in y
#' @param x an MID
#' @param y_carbons the number of carbons in the vector y
convolution_matrix <- function(x, y_carbons)
{
  x_carbons <- length(x) - 1
  # create empty matrix
  A <- matrix(0, x_carbons + y_carbons + 1, y_carbons + 1)
  # fill in columns
  for (i in 1:(y_carbons + 1)) {
    A[i:(i+x_carbons), i] <- x
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
  return(as.vector(convolution_matrix(x, length(y)-1) %*% y))
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

#
# 
# in cases of unequal carbons
#

#' Calculate similarity, possibly using convolution
#'
#' brings MIDs of peak x and peak y for experiment e from an MIData object (sim_data)
#' computes convolution in case of unequal C numbers
#' calculates similarity between MIDs
#' removes the first element that corresponds to M+0, if remove_m0 == T, before computing similarity.
#' @param sim_data the MIdata object
#' @param x peak index
#' @param y peak index
#' @param e experiment index
#' @param similarity similarity function
#' @param remove_m0 whether to exclude M+0 (the first element in MID vectors) from similarity calculation
#' @returns the convoluted MID vector, or NA if isotopoic enrichment is less than than tolerance
#' @export

conv_similarity <- function(sim_data, x, y, e, similarity, remove_m0 = F)
{
  # find in sim_data (an MIData object) the MID of metabolite with index x from experiment e
  mid_x <- get_avg_mid(sim_data, x, e)
  # number of carbon atoms of metabolite x
  n_atom_x <- length(mid_x) - 1
  
  # find in sim_data (an MIData object) the MID of metabolite with index y from experiment e
  mid_y <- get_avg_mid(sim_data, y, e)
  # number of carbon atoms of metabolite y
  n_atom_y <- length(mid_y) - 1
  
  if (n_atom_x == n_atom_y) {
    # in case of equal number, just calculate the similarity
    if(remove_m0 == TRUE)
    {
      mid_x <- mid_x[-1]
      mid_y <- mid_y[-1]
    }
    return(similarity(mid_x, mid_y))
  }
  
  else {
    # smaller and larger MIDs
    if(n_atom_x < n_atom_y) {
      s_mid <- mid_x
      l_mid <- mid_y
      carbon_diff <- n_atom_y - n_atom_x
    }
    else {
      s_mid <- mid_y
      l_mid <- mid_x
      carbon_diff <- n_atom_x - n_atom_y
    }
    
    # index of metabolites to convolute with
    conv_met_index <- get_peak_index_n_atoms(sim_data, carbon_diff)
    
    if (length(conv_met_index) > 0) {
      # get all possible convolutions
      convolutions <- lapply(
        conv_met_index,
        # for simulated data, get_avg_mid and get_mid would return the same MID from sim_data
        # for real data, this returns the average MID of all replicates
        function(m) convolute(get_avg_mid(sim_data, m, e), s_mid))
      
      # if sim_data contains C-13 corrected MIDs, remove M+0 from all convolutions
      if (remove_m0 == TRUE){
        convolutions <- lapply(convolutions, function(conv) return(conv[-1]))
        # remove M+0 also from the larger metabolite
        l_mid <- l_mid[-1]
      }
      
      # calculate similarities between the larger metabolite and all possible convolutions
      # without M+0
      similarities <- unlist(lapply(convolutions, similarity, l_mid))
      # assign zero similarity to missing values
      # missing values can occur after correcting unlabelled MIDs for natural 13C,
      # and the removal of M+0 results in zero vectors, thus zero division in similarity, e.g. cosine
      similarities[which(is.na(similarities))] <- 0
      # return the maximum similarity
      return(max(similarities))
    }
    else {
      # no matching metabolites to convolute with
      # return zero similarity
      return(0)
    }
  }
}
