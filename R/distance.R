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
#' experiment e from an MIData object (mi_data).
#' If x,y have unequal atom numbers, computes all convolutions x*z of the same
#' size as y (if y is the longer vector) and returns the largest similarity;
#' if similarity returns NA for some pair, the NA value is ignored when taking
#' maximum similarity; if there are no possible convolutions, returns 0
#'
#' @param mi_data the MIdata object
#' @param x peak index
#' @param y peak index
#' @param e experiment index
#' @param similarity A similarity function on MIDs.
#' @returns the resulting similarity value
#' @export

conv_similarity <- function(mi_data, x, y, e, similarity)
{
  # number of carbon atoms of metabolites x, y
  n_atom_x <- get_peak_n_atoms(mi_data, x)
  n_atom_y <- get_peak_n_atoms(mi_data, y)
  # ensure MID x is smaller or equal to MID y
  if(n_atom_x > n_atom_y) {
    return(conv_similarity(mi_data, y, x, e, similarity))
  }

  # find the MIDs of metabolite with index x, y from experiment e
  mid_x <- get_avg_mid(mi_data, x, e)
  mid_y <- get_avg_mid(mi_data, y, e)

  if (n_atom_x == n_atom_y) {
    # in case of equal number, just calculate the similarity
    return(similarity(mid_x, mid_y))
  }
  else {
    # x is strictly smaller than y
    # index of metabolites to convolute with
    carbon_diff <- n_atom_y - n_atom_x
    conv_met_index <- get_peak_index_n_atoms(mi_data, carbon_diff)

    if (length(conv_met_index) > 0) {
      # get all possible convolutions
      convolutions <- lapply(
        conv_met_index,
        function(m) convolute(get_avg_mid(mi_data, m, e), mid_x))

      # calculate similarities between the larger metabolite and all possible convolutions
      similarities <- unlist(lapply(convolutions, similarity, mid_y))
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


#' Calculate similarity matrix, based on a specific similarity/distance measure, for metabolites from an MIData object
#'
#' Calculates a similarity matrix based on the given similarity function for all MID pairs for
#' experiment e from an MIData object (midata).
#'
#' @param midata the MIdata object
#' @param e experiment index
#' @param similarity A similarity function on MIDs
#' @param remove_m0 whether to remove M+0 from MIDs (TRUE if MIDs are corrected for the natural 13C. FALSE by default)
#' @returns the resulting similarity matrix
#' @export

similarity_matrix <- function(midata, e, similarity, remove_m0 = FALSE)
{
  
  # c13 correct midata, and update similarity to remove M+0
  if (remove_m0 == TRUE)
  {
    midata <- midata_transform(midata, c13correct)
    similarity <- function(x,y) apply_no_m0(similarity, x,y)
  }
  
  # data dimensions
  n_metabolites <- length(midata$peak_ids)
  
  # met names
  met_names <- midata$peak_ids
  
  # create an empty matrix to be filled in with similarities
  sm <- matrix(NA, n_metabolites, n_metabolites)
  colnames(sm) <- midata$peak_ids
  rownames(sm) <- midata$peak_ids
  
  # loop over matrix elements (r2,d2) such that r2 <= d2
  for (x in 1:n_metabolites) 
  {
    for (y in x:n_metabolites) 
    {
      sm[x,y] <- sm[y,x] <- conv_similarity(midata, x, y, e, similarity)
    }
  }
  
  return(sm)
  
}
