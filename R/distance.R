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
  # compute norms
  x_norm <- sqrt(sum(x^2))
  y_norm <- sqrt(sum(y^2))
  # test for zero vector
  if(x_norm == 0 || y_norm == 0)
    return(NA)
  else
    return(sum(x * y) / (x_norm * y_norm))
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
  return(1 - cosine_sim(x,y))
}

#' Apply a function on MIDs x,y after removing M+0
#'
#' This simply removes the first component of each vector before applying f
#' @param f function to apply
#' @param x a vector
#' @param y a vector
#' @returns the result of calling f
#' @export
#'
apply_no_m0 <- function(f, x, y)
{
  return(f(x[-1], y[-1]))
}


calc_dot <- function(x, y) {1 - sum(x*y)}



#' Calculate similarity between average MIDs from an MIData object
#'
#' Calculates the given similarity function between MIDs for peak x and y for
#' experiment e from an MIData object (sim_data).
#' If x,y have unequal atom numbers, computes all convolutions x*z of the same
#' size as y (if y is the longer vector) and returns the largest similarity;
#' if similarity returns NA for some pair, the NA value is ignored when taking
#' maximum similarity.
#'
#' @param sim_data the MIdata object
#' @param x peak index
#' @param y peak index
#' @param e experiment index
#' @param similarity A similarity function on MIDs.
#' @returns the resulting similarity value
#' @export

conv_similarity <- function(sim_data, x, y, e, similarity)
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
        function(m) convolute(get_avg_mid(sim_data, m, e), s_mid))

      # calculate similarities between the larger metabolite and all possible convolutions
      similarities <- unlist(lapply(convolutions, similarity, l_mid))
      # remove any missing values
      similarities <- similarities[!is.na(similarities)]
      if(length(similarities) == 0) {
        # all NA, so return NA
        return(NA)
      }
      else {
        # return the maximum similarity
        return(max(similarities))
      }
    }
    else {
      # no matching metabolites to convolute with
      # return zero similarity
      return(0)
    }
  }
}
