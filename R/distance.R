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

#' Dot product similarity
#' @param x a vector
#' @param y a vector
#' @returns the dot (scalar) product x.y
#' @export
#'
dot_sim <- function(x, y)
{
  return(sum(x*y))
}

#' Dot product distance (not really a distance)
#' @param x a vector
#' @param y a vector
#' @returns a distance based on the dot (scalar) product x.y
#' This is always nonnegative if both vectors are MIDs, and zero if x = y
#' @export
#'
dot_dist <- function(x, y)
{
  return(sqrt(sum(x*x)*sum(y*y)) - sum(x*y))
}


#' The squared Euclidean distance (sum of squares)
#'
#' @param x a vector
#' @param y a vector
#' @returns the square of the Euclidean distance between x and y
#' @export
#'
euclidean_dist_sq <- function(x, y)
{
  diff <- x - y
  return(sum(diff*diff))
}

#' The Euclidean distance
#'
#' @param x a vector
#' @param y a vector
#' @returns the Euclidean distance between x and y
#' @export
#'
euclidean_dist <- function(x, y)
{
  return(sqrt(euclidean_dist_sq(x, y)))
}


#' Find maximum similarity between average MIDs from an MIData object
#'
#' Calculates the given similarity function between MIDs for peak x and y for
#' experiment e from an MIData object (mi_data).
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
  conv_reduce(mi_data, x, y, e, similarity, max_sim)
}

max_sim <- function(sim_values)
{
  if(length(sim_values) == 0)
    return(0)
  sim_no_na <- sim_values[!is.na(sim_values)]
  if(length(sim_no_na) > 0)
    return(max(sim_no_na))
  else
    return(NA)
}


#' Find minimum distance between average MIDs from an MIData object
#'
#' Calculates the given similarity function between MIDs for peak x and y for
#' experiment e from an MIData object (mi_data).
#'
#' @param mi_data the MIdata object
#' @param x peak index
#' @param y peak index
#' @param e experiment index
#' @param distance A similarity function on MIDs.
#' @param default Distance value to use when there is no possible convolution
#' @returns the resulting similarity value
#' @export
#'
conv_distance <- function(mi_data, x, y, e, distance, default = Inf)
{
  conv_reduce(mi_data, x, y, e, distance,
              function(dist_values) min_dist(dist_values, default))
}

min_dist <- function(dist_values, default = Inf)
{
  if(length(dist_values) == 0)
    return(default)
  return(min(dist_values, na.rm = TRUE))
}

# could perhaps be moved to midata.R
get_avg_mids_by_size <- function(mi_data, n_atoms, e)
{
  index <- get_peak_index_n_atoms(mi_data, n_atoms)
  return(sapply(index, function(i) get_avg_mid(mi_data, i, e)))
}


#' Calculates f(x, y) between MIDs for peak x and y for
#' experiment e from an MIData object (mi_data)
#' If x,y have unequal atom numbers, computes f(x*z, y) for all atoms of the same
#' size as y (if y is the longer vector) and returns g( ... ) of the list of
#' value of f. For example f = max gives the maximum value. The function g
#' must handle empty lists g(c()) in case there are no possible convolutions.
#'
#' @param mi_data the MIdata object
#' @param x peak index
#' @param y peak index
#' @param e experiment index
#' @param f A function f(x, y) taking two MIDs.
#' @param g a function g taking a vector of values f1, f2, ...
#' @returns the resulting value g(f1, f2, ...)
#' @export

conv_reduce <- function(mi_data, x, y, e, f, g)
{
  # number of carbon atoms of metabolites x, y
  n_atom_x <- get_peak_n_atoms(mi_data, x)
  n_atom_y <- get_peak_n_atoms(mi_data, y)
  # ensure MID x is smaller or equal to MID y
  if(n_atom_x > n_atom_y) {
    return(conv_reduce(mi_data, y, x, e, f, g))
  }
  # find the MIDs of metabolite with index x, y from experiment e
  mid_x <- get_avg_mid(mi_data, x, e)
  mid_y <- get_avg_mid(mi_data, y, e)

  if (n_atom_x == n_atom_y) {
    # equal numbers just calculate f
    return(g(c(f(mid_x, mid_y))))
  }
  else {
    # x is strictly smaller than y
    # get MIDs of metabolites z to convolute with
    n_atom_z <- n_atom_y - n_atom_x
    mids_z <- get_avg_mids_by_size(mi_data, n_atom_z, e)
    if(length(mids_z) > 0) {
      # compute all convolutions x*z for each z
      convolutions <- convolute_cols(mid_x, mids_z)
      # calculate f between the larger metabolite and all possible convolutions
      f_values <- apply(convolutions, MARGIN = 2, f, mid_y)
      # return the function g
      return(g(f_values))
    }
    else {
      # no matching metabolites to convolute with
      return(g(c()))
    }
  }
}


#' Same as conv_reduce, but for all pairs x*y, returns a matrix.
#' this avoids recalculating convolutions
#' @param mi_data the MIdata object
#' @param e experiment index
#' @param f A function f(x, y) taking two MIDs.
#' @param g a function g taking a vector of values f1, f2, ...
#' @returns the matrix of g(f(x,y) ...) values for all x,y
#' @export
conv_reduce_all <- function(mi_data, e, f, g)
{
  # allocate square matrix
  n_met <- length(mi_data$peak_ids)
  result <- matrix(0, nrow = n_met, ncol = n_met)
  # unique metabolite sizes
  n_atoms <- unique(mi_data$peak_n_atoms)
  n_atoms <- n_atoms[order(n_atoms)]
  # loop over x
  for(i in 1:length(n_atoms)) {
    n_atoms_x <- n_atoms[[i]]
    for(j in i : length(n_atoms)) {
      n_atoms_y <- n_atoms[[j]]
      # get all MIDs y
      y_index <- get_peak_index_n_atoms(mi_data, n_atoms_y)
      mids_y <- sapply(y_index, function(i) get_avg_mid(mi_data, i, e))

      # metabolites z to convolute with x
      n_atoms_z <- n_atoms_y - n_atoms_x
      mids_z <- get_avg_mids_by_size(mi_data, n_atoms_z, e)

      for(x in get_peak_index_n_atoms(mi_data, n_atoms_x)) {
        # MIDs of metabolite with index x
        mid_x <- get_avg_mid(mi_data, x, e)

        if(n_atoms_x == n_atoms_y) {
          # for metabolites y of same size as x, just calculate f(x,y)
          result[x, y_index] <- apply(mids_y, MARGIN = 2,
                                      function(y) g(c(f(mid_x, y))))
        }
        else {
          if(length(mids_z) > 0) {
            # compute all convolutions x*z for each z
            convolutions <- convolute_cols(mid_x, mids_z)
            # for each y, calculate g(f(x*z, y)) over all convolutions x*z
            result[x, y_index] <- apply(mids_y, MARGIN = 2,
                     function(y) g(apply(convolutions, MARGIN = 2, f, y)))
          }
          else {
            # no metabolites z to convolute x with
            result[x, y_index] <- rep(g(c()), n = length(y_index))
          }
        }
        result[y_index, x] <- result[x, y_index]
      }
    }
  }
  return(result)
}




