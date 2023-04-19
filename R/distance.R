
# distance and similarity measures
#

#' Apply a function on MIDs x,y after removing M+0
#'
#' This simply removes the first component of each vector before applying f
#' @param f function to apply
#' @param x a vector
#' @param y a vector
#' @returns the result of calling f
#' @export
#'
apply_no_m0 <- function(f, x, y) {
  return(f(x[-1], y[-1]))
}

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
#' This is always nonnegative if both vectors are MIDs, and zero if x = y
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
#' This is always nonnegative if both vectors are MIDs, and zero if x = y
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
  sum(abs(x - y))
}

#' @export
canberra_distance <- function(x, y) {
  sum(abs(x - y) / (abs(x) + abs(y)))
}

#' @export
bray_curtis_distance <- function(x, y) {
  sum(abs(x - y)) / sum(abs(x + y))
}


# function adding z_index vector as attribute to f values
f_attr <- function(mid_y, mids_xz, z_index, f)
{
  f_values <- apply(mids_xz, MARGIN = 2, f, mid_y)
  attr(f_values, "index") <- z_index
  return(f_values)
}


# alternative, use a list with $value and $index instead of attribute
f_list <- function(mid_y, mids_xz, z_index, f)
{
  f_values <- apply(mids_xz, MARGIN = 2, f, mid_y)
  return(list(values = f_values, index = z_index))
}


#' Calculates a convolution-based measure between MIDs for peak x and y for
#' experiment e from an MIData object (mi_data).
#'
#' The function f(x, y) must be symmetric in x, y and return a real value
#' or NA for any x, y with the same number of atoms.
#' If x, y have the same number of atoms, f(x,y) is returned.
#'
#' If x is smaller than y, this function computes f(x*z, y) for all
#' peaks z such that x*z is the same size as y (and vice versa), and calls g()
#' on the resulting vector of values from f.
#' For example, setting f = max gives the maximum value.
#' The function g must accept an empty lists g(c()) in case there are no
#' possible convolutions, and handle all values produced by f.
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
  if (n_atom_x > n_atom_y) {
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
    z_index <- get_peak_index_n_atoms(mi_data, n_atom_z)
    mids_z <- sapply(z_index, function(i) get_avg_mid(mi_data, i, e))
    if (length(mids_z) > 0) {
      # compute all convolutions x*z for each z
      mids_xz <- convolute_cols(mid_x, mids_z)
      # calculate f between y and all x*y and add indices
      f_values <- f_attr(mid_y, mids_xz, z_index, f)
      # return the function g
      return(g(f_values))
    } else {
      # no matching metabolites to convolute with
      return(g(c()))
    }
  }
}


#' Same as conv_reduce, but calculates f(x, y) pairwise for all peaks x and y
#' from an MIData object (mi_data)
#'
#' The function g must handle empty lists g(c()) in case there are no possible convolutions.
#' In that case, g_select will not be applied.
#'
#' @param mi_data the MIdata object
#' @param e experiment index
#' @param f A function f(x, y) taking two MIDs.
#' @param g a function g taking a vector of values f1, f2, ...
#' @returns the matrix of g(f(x,y) ...) values for all x,y
#' @export
conv_reduce_all <- function(mi_data, e, f, g)
{
  n_met <- length(mi_data$peak_ids)
  # allocate matrices
  conv_values <- matrix(as.double(NA), n_met, n_met)
  conv_index <- matrix(as.integer(NA), n_met, n_met)

  # unique metabolite sizes, sorted
  n_atoms <- sort(unique(mi_data$peak_n_atoms))
  # for each metabolite size for x
  for(i in 1:length(n_atoms)) {
    n_atoms_x <- n_atoms[[i]]
    # get all MIDs for x
    x_index <- get_peak_index_n_atoms(mi_data, n_atoms_x)

    # for each metabolites size for y, larger than x
    for(j in i:length(n_atoms)) {
      n_atoms_y <- n_atoms[[j]]
      # get all MIDs y of this size
      y_index <- get_peak_index_n_atoms(mi_data, n_atoms_y)
      # compute this block
      if(n_atoms_x == n_atoms_y) {
        block <- conv_reduce_block_equal(mi_data, e, f, g, x_index, y_index)
      }
      else {
        n_atoms_z = n_atoms_y - n_atoms_x
        z_index <- get_peak_index_n_atoms(mi_data, n_atoms_z)
        block <- conv_reduce_block(mi_data, e, f, g, x_index, y_index, z_index)
      }
      stopifnot(is.matrix(block$values))
      stopifnot(is.matrix(block$index))
      # copy block to full matrices
      conv_values[x_index, y_index] <- block$values
      conv_values[y_index, x_index] <- t(block$values)
      conv_index[x_index, y_index] <- block$index
      conv_index[y_index, x_index] <- t(block$index)
    }
  }
  return(with_attr(conv_values, "index", conv_index))
}



# apply a selection function g and return a list of
# g values and corresponding indices
g_list <- function(mids_y, mids_xz, z_index, f, g)
{
  n_y <- dim(mids_y)[2]
  g_values <- rep(as.double(NA), n_y)
  g_index <- rep(as.integer(NA), n_y)
  for(i in 1:n_y) {
    g_result <- g(f_attr(mids_y[, i], mids_xz, z_index, f))
    # store values and index separately
    g_values[i] <- without_attr(g_result)
    g_index[i] <- attr(g_result, "index")
  }
  # add index vector as attribute
  return(list(values = g_values, index = g_index))
}



#' compute conv_reduce for one matrix block where all x are the same size
#' and all y are the same size, but y is larger than x
#'
#' @param mids_x matrix of MIDs for metabolites x, each column one MID
#' @param mids_y matrix of MIDs for metabolites y
#' @param mids_z matrix of MIDs for metabolites z, may be empty
#'
conv_reduce_block <- function(mi_data, e, f, g, x_index, y_index, z_index)
{
  n_x <- length(x_index)
  n_y <- length(y_index)
  ## allocate matrices
  block = list(values = matrix(as.double(NA), n_x, n_y),
               index = matrix(as.integer(NA), n_x, n_y))
  if(length(z_index) > 0) {
    # get MID matrices
    mids_x <- sapply(x_index, function(i) get_avg_mid(mi_data, i, e))
    stopifnot(is.matrix(mids_x))
    mids_y <- sapply(y_index, function(i) get_avg_mid(mi_data, i, e))
    stopifnot(is.matrix(mids_y))
    mids_z <- sapply(z_index, function(i) get_avg_mid(mi_data, i, e))
    stopifnot(is.matrix(mids_z))

    for(i in 1:n_x) {
      # compute all convolutions x*z for each z
      mids_xz <- convolute_cols(mids_x[, i], mids_z)
      stopifnot(is.matrix(mids_xz))
      # calculate g(f(x*z, y) ...) for all y (rows) and all convolutions x*z
      g_result <- g_list(mids_y, mids_xz, z_index, f, g)
      # NOTE: is.vector(x) is FALSE when x is a vector but has attributes!
      stopifnot(is.vector(g_result$values))
      # store value and attributes separately
      block$values[i, ] <- g_result$values
      block$index[i, ] <- g_result$index
    }
  }
  # else no metabolites z to convolute x with
  return(block)
}


# special case of the above, for a block where all x AND y are the same size
# in this case there is nothing to convolute with
conv_reduce_block_equal <- function(mi_data, e, f, g, x_index, y_index)
{
  n_x <- length(x_index)
  n_y <- length(y_index)
  # allocate matrices
  block = list(values = matrix(as.double(NA), n_x, n_y),
               index = matrix(as.integer(NA), n_x, n_y))

  mids_x <- sapply(x_index, function(i) get_avg_mid(mi_data, i, e))
  stopifnot(is.matrix(mids_x))
  mids_y <- sapply(y_index, function(i) get_avg_mid(mi_data, i, e))
  stopifnot(is.matrix(mids_y))

  # calculate g(f(x,y)) for each x,y (no attributes in this case)
  for(i in 1:n_x) {
    block$values[i, ] <- apply(mids_y, MARGIN = 2,
                               function(mid_y) g(c(f(mids_x[, i], mid_y))))
  }
  return(block)
}

conv_reduce_block_equal_all_experiments <- function(mi_data, f, g, x_index, y_index)
{
  n_x <- length(x_index)
  n_y <- length(y_index)
  # allocate matrices
  block = list(values = matrix(as.double(NA), n_x, n_y),
               index = matrix(as.integer(NA), n_x, n_y))
  
  mids_x <- lapply(1:length(mi_data$experiments), 
                   function(e) sapply(x_index, function(i) get_avg_mid(mi_data, i, e)))
  stopifnot(is.matrix(mids_x[[1]]))
  mids_y <- lapply(1:length(mi_data$experiments), 
                   function(e) sapply(y_index, function(i) get_avg_mid(mi_data, i, e)))
  stopifnot(is.matrix(mids_y[[1]]))
  
  # calculate g(f(x,y)) for each x,y (no attributes in this case)
  for(i in 1:n_x) {
    block$values[i, ] <- apply(mids_y, MARGIN = 2,
                               function(mid_y) g(c(f(mids_x[, i], mid_y))))
  }
  return(block)
}

euclidean_dist_sq_from_list <- function(x, y)
{
  # ensure equal number of experiments
  stopifnot(length(x) == length(y))
  
  for (i in 1:length(x)){
    
  }
  stopifnot(is.vector(x))
  stopifnot(is.vector(y))
  stopifnot(length(x) == length(y))
  diff <- x - y
  return(sum(diff * diff))
}

#' @export
conv_reduce_all_test <- function(mi_data, e, f, g, get_middle_met_matrix = F, g_select, z_threshold = 0.0107) {
  # allocate square matrix
  n_met <- length(mi_data$peak_ids)
  result <- matrix(0, nrow = n_met, ncol = n_met)
  if (get_middle_met_matrix == T) {
    middle_met_matrix <- matrix(NA, nrow = n_met, ncol = n_met)
  }

  # unique metabolite sizes
  n_atoms <- unique(mi_data$peak_n_atoms)
  n_atoms <- n_atoms[order(n_atoms)]
  # loop over x
  for (i in 1:length(n_atoms)) {
    n_atoms_x <- n_atoms[[i]]
    for (j in i:length(n_atoms)) {
      n_atoms_y <- n_atoms[[j]]
      # get all MIDs y
      y_index <- get_peak_index_n_atoms(mi_data, n_atoms_y)
      mids_y <- sapply(y_index, function(i) get_avg_mid(mi_data, i, e))

      # metabolites z to convolute with x
      n_atoms_z <- n_atoms_y - n_atoms_x
      mids_z <- as.matrix(get_avg_mids_by_size(mi_data, n_atoms_z, e))


      for (x in get_peak_index_n_atoms(mi_data, n_atoms_x)) {
        # MIDs of metabolite with index x
        mid_x <- get_avg_mid(mi_data, x, e)

        if (n_atoms_x == n_atoms_y) {
          # for metabolites y of same size as x, just calculate f(x,y)
          result[x, y_index] <- result[y_index, x] <- apply(mids_y,
                                                            MARGIN = 2,
                                                            function(y) g(c(f(mid_x, y)))
          )
        } else {
          if (length(mids_z) > 0){
            # compute isotopic enrichment
            z_enrichment <- unlist(apply(mids_z, 2, isotopic_enrichment))

            if (length(which(z_enrichment > z_threshold)) > 0) {
              # compute all convolutions x*z for each z
              convolutions <- convolute_cols(mid_x, as.matrix(mids_z[, which(z_enrichment > z_threshold)]))
              # get indices of all 13-C enriched z
              z_ind <- get_peak_index_n_atoms(mi_data, n_atoms_z)[which(z_enrichment > z_threshold)]

              # for each y, calculate g(f(x*z, y)) over all convolutions x*z
              # get pairwise measures for each convolution
              pm_conv <- apply(mids_y, 2, function(y) apply(convolutions, MARGIN = 2, f, y))

              if (length(z_ind) == 1) {
                result[x, y_index] <- result[y_index, x] <- pm_conv
                if (get_middle_met_matrix == T) {
                  middle_met_matrix[x, y_index] <- middle_met_matrix[y_index, x] <- z_ind
                }
              } else {
                result[x, y_index] <- result[y_index, x] <- apply(pm_conv, 2, g)
                if (get_middle_met_matrix == T) {
                  middle_met_matrix[x, y_index] <- middle_met_matrix[y_index, x] <- z_ind[apply(pm_conv, 2, g_select)]
                }
              }
            }
          }
          else {
            # no metabolites z to convolute x with
            result[x, y_index] <- result[y_index, x] <- rep(suppressWarnings(g(c())), n = length(y_index))
            if (get_middle_met_matrix == T) {
              middle_met_matrix[x, y_index] <- middle_met_matrix[y_index, x] <- NA
            }
          }
        }
      }
    }
  }
  if (get_middle_met_matrix == T) {
    return(list(result, middle_met_matrix))
  } else {
    return(result)
  }
}


#' @export
conv_reduce_all_baseline <- function(mi_data, e, f, binom_prob = 0.0107) {
  # allocate square matrix
  n_met <- length(mi_data$peak_ids)
  result <- matrix(0, nrow = n_met, ncol = n_met)

  # unique metabolite sizes
  n_atoms <- unique(mi_data$peak_n_atoms)
  n_atoms <- n_atoms[order(n_atoms)]

  # loop over x
  for (i in 1:length(n_atoms)) {
    n_atoms_x <- n_atoms[[i]]
    for (j in i:length(n_atoms)) {
      n_atoms_y <- n_atoms[[j]]
      # get all MIDs y
      y_index <- get_peak_index_n_atoms(mi_data, n_atoms_y)
      mids_y <- sapply(y_index, function(i) get_avg_mid(mi_data, i, e))

      # metabolites z to convolute with x
      n_atoms_z <- n_atoms_y - n_atoms_x

      for (x in get_peak_index_n_atoms(mi_data, n_atoms_x)) {
        # MIDs of metabolite with index x
        mid_x <- get_avg_mid(mi_data, x, e)

        if (n_atoms_x == n_atoms_y) {
          # for metabolites y of same size as x, just calculate f(x,y)
          result[x, y_index] <- result[y_index, x] <- apply(mids_y,
                                                            MARGIN = 2,
                                                            function(y) f(mid_x, y)
          )
        } else {
          # binomial distribution
          bino <- binomvals(0:n_atoms_z, n_atoms_z, binom_prob)
          # convolute with unlabelled mid
          convolution <- convolute(mid_x, bino)
          # for each y, calculate g(f(x*z, y))
          result[x, y_index] <- result[y_index, x] <- apply(mids_y, 2, f, convolution)
        }
      }
    }
  }
  return(result)
}



#' Calculate similarity matrix, based on parameters specified in an InputData object
#'
#' Calculates a similarity matrix based on the given similarity function for all peak pairs
#' for experiment e from an MIData object (midata).
#' If a list of reactions is input to input_data$reaction_data, the given similarity function will only be applied to peak pairs
#' whose formula or mass difference (depending on "input_data$reaction_restriction") can be matched an allowed reaction
#'
#' @param e experiment index that corresponds to the same index in MIData
#' @param input_data InputData object (for more detail see ...)
#' @export

pairwise_matrix_all_ds <- function(e, input_data) {
  experiment <- input_data$midata$experiments[e]

  # print progress
  cat("Computing ", input_data$measure, " matrix for experiment: ", experiment, "\n")
  if (input_data$reaction_restriction == F) {
    cat("for all possible convolutions in data", "\n")
  } else if (input_data$reaction_restriction == "mass") {
    cat("for a restriction on mass difference", "\n")
  } else {
    cat("for a restriction on formula difference", "\n")
  }

  # data dimensions
  n_metabolites <- length(input_data$midata$peak_ids)
  # met names
  met_names <- input_data$midata$peak_ids

  # compute the pairwise matrix
  pairwise_matrix <- conv_reduce_all(
    input_data$midata, e,
    input_data$fun,
    input_data$perfection,
    eval(parse(text = input_data$get_middle_met_matrix)),
    input_data$g_select
  )


  # assigning column names
  if (input_data$get_middle_met_matrix == T) {
    colnames(pairwise_matrix[[1]]) <- rownames(pairwise_matrix[[1]]) <-
      colnames(pairwise_matrix[[2]]) <- rownames(pairwise_matrix[[2]]) <-
      input_data$midata$peak_ids
  } else {
    colnames(pairwise_matrix) <- rownames(pairwise_matrix) <- input_data$midata$peak_ids
  }

  # return the final matrix (or matrices)
  return(pairwise_matrix)
}

#' @export
pairwise_matrix_all_ds_test <- function(e, input_data) {
  experiment <- input_data$midata$experiments[e]

  # print progress
  cat("Computing ", input_data$measure, " matrix for experiment: ", experiment, "\n")
  if (input_data$reaction_restriction == F) {
    cat("for all possible convolutions in data", "\n")
  } else if (input_data$reaction_restriction == "mass") {
    cat("for a restriction on mass difference", "\n")
  } else {
    cat("for a restriction on formula difference", "\n")
  }


  # data dimensions
  n_metabolites <- length(input_data$midata$peak_ids)
  # met names
  met_names <- input_data$midata$peak_ids

  # compute the pairwise matrix
  pairwise_matrix <- conv_reduce_all_test(
    input_data$midata, e,
    input_data$fun,
    input_data$perfection,
    eval(parse(text = input_data$get_middle_met_matrix)),
    input_data$g_select
  )

  # assigning column names
  if (input_data$get_middle_met_matrix == T) {
    colnames(pairwise_matrix[[1]]) <- rownames(pairwise_matrix[[1]]) <-
      colnames(pairwise_matrix[[2]]) <- rownames(pairwise_matrix[[2]]) <-
      input_data$midata$peak_ids
  } else {
    colnames(pairwise_matrix) <- rownames(pairwise_matrix) <- input_data$midata$peak_ids
  }

  # return the final matrix (or matrices)
  return(pairwise_matrix)
}

#' @export
pairwise_matrix_all_ds_baseline <- function(e, input_data) {
  experiment <- input_data$midata$experiments[e]

  # print progress
  cat("Computing ", input_data$measure, " matrix for experiment: ", experiment, "\n")
  if (input_data$reaction_restriction == F) {
    cat("for all possible convolutions in data", "\n")
  } else if (input_data$reaction_restriction == "mass") {
    cat("for a restriction on mass difference", "\n")
  } else {
    cat("for a restriction on formula difference", "\n")
  }


  # data dimensions
  n_metabolites <- length(input_data$midata$peak_ids)
  # met names
  met_names <- input_data$midata$peak_ids

  # compute the pairwise matrix
  pairwise_matrix <- conv_reduce_all_baseline(
    input_data$midata, e,
    input_data$fun
  )

  # assigning column names
  colnames(pairwise_matrix) <- rownames(pairwise_matrix) <- input_data$midata$peak_ids

  # return the final matrix (or matrices)
  return(pairwise_matrix)
}

#' Weighs the pairwise measure by the enrichment scores of metabolites involved.
#' Assumes distances only.
#' In case of equal C, the pairwise distance is divided by isotopic_enrichment(mid_x)*isotopic_enrichment(mid_y) for pair (x,y)
#' In case of unequal C, it is divided by isotopic_enrichment(mid_x)*isotopic_enrichment(mid_z)*isotopic_enrichment(mid_y) for pair (x*z,y)
#' Returns a penalized pairwise distance matrix.
#'
#' @param pairwise_matrix distance matrix
#' @param middle_met_matrix matrix of middle metabolites that were chosen as the best convolutions
#' @param experiment_matrix matrix of the same size, whose cells denote the experiment that was the output of g_select()
#' @param mi_data the MIData object
#' @export
weight_pm_by_enrichment <- function(pairwise_matrix, middle_met_matrix, experiment_matrix, mi_data) {
  # create an empty matrix to be filled in by weighted pairwise measures
  weighted_pm <- matrix(NA, nrow(pairwise_matrix), ncol(pairwise_matrix))
  colnames(weighted_pm) <- rownames(weighted_pm) <- colnames(pairwise_matrix)
  # loop over rows and columns
  for (i in 1:nrow(pairwise_matrix)) {
    for (j in i:nrow(pairwise_matrix)) {
      # pairwise measure should be non infinite
      if (is.infinite(pairwise_matrix[i, j]) == F) {
        ie_1 <- isotopic_enrichment(get_avg_mid(mi_data, i, experiment_matrix[i, j]))
        ie_2 <- isotopic_enrichment(get_avg_mid(mi_data, j, experiment_matrix[i, j]))
        # for equal C pairs multiply enrichments of the two MIDs
        if (is.na(middle_met_matrix[i, j]) == T) {
          weighted_pm[i, j] <- weighted_pm[j, i] <- pairwise_matrix[i, j] / (ie_1 * ie_2)
        } else {
          # for different C number we need the middle metabolite
          ie_3 <- isotopic_enrichment(get_avg_mid(mi_data, middle_met_matrix[i, j], experiment_matrix[i, j]))
          weighted_pm[i, j] <- weighted_pm[j, i] <- pairwise_matrix[i, j] / (ie_1 * ie_2 * ie_3)
        }
      }
    }
  }
  return(weighted_pm)
}


#' Filters the pairwise distance matrix by only including the top percent distances from the CDF.
#' Filtering is performed for each row.
#' Assumes distances only.
#'
#' @param pairwise_matrix a pairwise distance matrix
#' @param percentile top X percent of the predictions will be kept where X is 1 per default.
#' @export
filter_pairwise_matrix <- function(pairwise_matrix, percentile = 0.01)
{
  # create an empty matrix to be filled in only by those who pass the filtering criteria
  filtered_pm <- matrix(NA, nrow(pairwise_matrix), ncol(pairwise_matrix))
  colnames(filtered_pm) <- rownames(filtered_pm) <- colnames(pairwise_matrix)

  # make selections for each row
  for (i in 1:nrow(pairwise_matrix))
  {
    # get infinite and NA indices
    non_inf_ind <- which(is.infinite(pairwise_matrix[i, ]) == F & is.na(pairwise_matrix[i, ]) == F)
    # exclude infinites and NAs
    non_inf_vec <- pairwise_matrix[i, ][non_inf_ind]

    # compute the percentile (top 1%)
    threshold <- as.numeric(quantile(non_inf_vec, probs = percentile))

    # fill in
    filtered_pm[i, as.numeric(non_inf_ind[which(non_inf_vec <= threshold)])] <- pairwise_matrix[i, as.numeric(non_inf_ind[which(non_inf_vec <= threshold)])]
  }
  return(filtered_pm)
}


#' @export
filter_pairwise_matrix_global <- function(pairwise_matrix, percentile = 0.01) {
  # remove diagonals
  diag(pairwise_matrix) <- NA

  # create an empty matrix to be filled in only by those who pass the filtering criteria
  filtered_pm <- matrix(NA, nrow(pairwise_matrix), ncol(pairwise_matrix))
  colnames(filtered_pm) <- rownames(filtered_pm) <- colnames(pairwise_matrix)

  # get infinite and NA indices
  non_inf_ind <- which(is.infinite(pairwise_matrix) == F & is.na(pairwise_matrix) == F)
  # exclude infinites and NAs
  non_inf_vec <- pairwise_matrix[non_inf_ind]
  # compute the percentile (top 1%)
  threshold <- as.numeric(quantile(non_inf_vec, probs = percentile))
  # fill in
  filtered_pm[as.numeric(non_inf_ind[which(non_inf_vec <= threshold)])] <- pairwise_matrix[as.numeric(non_inf_ind[which(non_inf_vec <= threshold)])]

  return(filtered_pm)
}


#' @export
subsample_pairwise_matrix <- function(pairwise_matrix, percentile = 0.01)
{
  # create an empty matrix
  subsampled_pm <- matrix(NA, nrow(pairwise_matrix), ncol(pairwise_matrix))
  colnames(subsampled_pm) <- rownames(subsampled_pm) <- colnames(pairwise_matrix)

  # make selections for each row
  for (i in 1:nrow(pairwise_matrix))
  {
    # get infinite and NA indices
    non_inf_ind <- which(is.infinite(pairwise_matrix[i, ]) == F & is.na(pairwise_matrix[i, ]) == F)

    # compute the percentile (top 1%)
    ind <- sample(non_inf_ind, round(length(non_inf_ind)*percentiles[p]*0.01))

    # fill in
    subsampled_pm[i, ind] <- pairwise_matrix[i, ind]
  }

  return(subsampled_pm)
}


#' Combines pairwise matrices across multiple experiments by their squared sum.
#' Note that this cannot be used if the middle metabolite matrix and experiment index matrix are to be tracked.
#'
#' @param pairwise_matrices a list of pairwise matrices across different experiments
#' @export
combine_sqrt_sum <- function(pairwise_matrices)
{
  return(sqrt(Reduce("+", pairwise_matrices)))
}


#' Combines pairwise matrices across multiple experiments by their sum.
#' Note that this cannot be used if the middle metabolite matrix and experiment index matrix are to be tracked.
#'
#' @param pairwise_matrices a list of pairwise matrices across different experiments
#' @export
combine_sum <- function(pairwise_matrices)
{
  return(Reduce("+", pairwise_matrices))
}


#' Returns the metabolite index which was chosen by fun.
#' It is supposed to be used in the combine() function.
#' @param vector a row from a pairwise matrix
#' @param fun a function like max() or min()
#' @export
get_fun_index <- function(vector, fun)
{
  vector <- as.vector(vector)
  non_na_inf_ind <- which(is.na(vector) == F & is.infinite(vector) == F)
  new_vector <- vector[non_na_inf_ind]

  return(non_na_inf_ind[which(new_vector == fun(new_vector))][1])
}


#' Combines pairwise matrices across multiple experiments based on the fun function.
#' Suitable for cases where the middle metabolite matrix is known and experiment index
#' matrix is to be tracked.
#'
#' @param pairwise_matrices a list of pairwise matrices across different experiments
#' @param middle_met_matrices a list of middle metabolite matrices across different
#' experiments, where a middle metabolite denotes the metabolite that was chosen
#' as the best convolution for a given metabolite pair of unequal carbon number
#'
#' @export
combine <- function(pairwise_matrices, middle_met_matrices, fun)
{
  pairwise_vec <- do.call(rbind.data.frame, lapply(pairwise_matrices, as.vector))
  middle_met_vec <- do.call(rbind.data.frame, lapply(middle_met_matrices, as.vector))

  index <- as.vector(unlist(apply(pairwise_vec, 2, get_fun_index, fun)))

  pm <- matrix(
    unlist(lapply(1:length(index),
                  function(x, pairwise_vec, index) pairwise_vec[[x]][index[[x]][1]],
                  pairwise_vec, index)),
    ncol = ncol(pairwise_matrices[[1]]),
    byrow = F
  )

  mmm <- matrix(
    unlist(lapply(1:length(index),
                  function(x, pairwise_vec, index) middle_met_vec[[x]][index[[x]][1]],
                  vector, index)),
    ncol = ncol(pairwise_matrices[[1]]),
    byrow = F
  )

  ei <- matrix(index,
    ncol = ncol(pairwise_matrices[[1]]),
    byrow = F
  )

  return(list(pm, mmm, ei))
}


#' @export
combine_test <- function(pairwise_matrices, middle_met_matrices, fun, input) {

  pairwise_vec <- do.call(rbind.data.frame, lapply(pairwise_matrices, as.vector))
  middle_met_vec <- do.call(rbind.data.frame, lapply(middle_met_matrices, as.vector))

  index <- as.vector(unlist(apply(pairwise_vec, 2, get_fun_index, fun)))

  pm <- matrix(unlist(lapply(1:length(index), function(x, pairwise_vec, index) pairwise_vec[[x]][index[[x]][1]], pairwise_vec, index)),
               ncol = ncol(pairwise_matrices[[1]]),
               byrow = F
  )

  mmm <- matrix(unlist(lapply(1:length(index), function(x, pairwise_vec, index) middle_met_vec[[x]][index[[x]][1]], vector, index)),
                ncol = ncol(pairwise_matrices[[1]]),
                byrow = F
  )

  ei <- matrix(index,
               ncol = ncol(pairwise_matrices[[1]]),
               byrow = F
  )
  return(list(pm, mmm, ei))
}

#' @export
convert_to_edge_list <- function(pairwise_matrix, middle_met_matrix, input, percentile){

  # filter pairwise matrix for the given percentile
  pairwise_matrix <- filter_pairwise_matrix_global(pairwise_matrix, percentile*0.01)
  
  edge_list <- data.frame(metabolite_1 = NA, met_1_attribute = NA,
                          metabolite_2 = NA, met_2_attribute = NA,
                          distance = NA, distance_attribute = NA)
  for (r2 in 1:nrow(pairwise_matrix)){
    for (d2 in r2:nrow(pairwise_matrix)){
  
      if (is.na(pairwise_matrix[r2, d2]) == F)
      {
        peak_1 <- input$midata$peak_ids[r2]
        nr_carbon_1 <- input$midata$peak_n_atoms[r2]
        peak_2 <- input$midata$peak_ids[d2]
        nr_carbon_2 <- input$midata$peak_n_atoms[d2]
  
        # in case of equal carbons
        if (nr_carbon_1 == nr_carbon_2)
        {
          edge_list <- rbind(edge_list,
                             data.frame(metabolite_1 = peak_1, met_1_attribute = "metabolite",
                                        metabolite_2 = peak_2, met_2_attribute = "metabolite",
                                        distance = pairwise_matrix[r2, d2],
                                        distance_attribute = "equal_carbon"))
        }
  
        # if metabolite 1 has less carbons than metabolite 2
        # peak_1 <- peak_convolution & peak_middle <- peak_convolution & peak_convolution <- peak_2
        else if (nr_carbon_1 < nr_carbon_2)
        {
          peak_middle <- input$midata$peak_ids[middle_met_matrix[r2,d2]]
          peak_convolution <- paste(r2, d2, sep = "_")
          edge_list <- rbind(edge_list,
                             data.frame(metabolite_1 = c(peak_1, peak_middle, peak_convolution),
                                        met_1_attribute = c("metabolite", "metabolite", "convolution"),
                                        metabolite_2 = c(peak_convolution, peak_convolution, peak_2),
                                        met_2_attribute = c("convolution", "convolution", "metabolite"),
                                        distance = pairwise_matrix[r2, d2], distance_attribute = "convolution"))
        }
  
        # if metabolite 1 has more carbons than metabolite 2
        # peak_2 <- peak_convolution & peak_middle <- peak_convolution & peak_convolution <- peak_1
        else
        {
          peak_middle <- input$midata$peak_ids[middle_met_matrix[r2,d2]]
          peak_convolution <- paste(r2, d2, sep = "_")
          edge_list <- rbind(edge_list,
                             data.frame(metabolite_1 = c(peak_2, peak_middle, peak_convolution),
                                        met_1_attribute = c("metabolite", "metabolite", "convolution"),
                                        metabolite_2 = c(peak_convolution, peak_convolution, peak_1),
                                        met_2_attribute = c("convolution", "convolution", "metabolite"),
                                        distance = pairwise_matrix[r2, d2], distance_attribute = "convolution"))
        }
      }
  
  
    }
  }
  # remove the first NA row
  edge_list <- edge_list[-1,]
  
  return(edge_list)
}

#' @export
# compute distances for equal carbon metabolites, and make a selection based on g_select
calc_distance_equal_c <- function(mid_x, mid_y, f){
  
  cd <- c()
  for (e in 1:length(mid_x)){
    cd <- c(cd, f(mid_x[[e]], mid_y[[e]]))
  }
  # returns 2 elements: (1) the final distance, (2) the middle metabolite which was chosen by g_select
  return(list(sum(cd), NA))
}

#' @export
# compute distances for unequal carbon metabolites, and make a selection based on g_select
calc_distance_unequal_c <- function(mid_x, mid_y, mid_z, z_index, f, g_select){
  sums_cd <- c()
  for (z in 1:length(mid_z)){
    cd <- c()
    for (e in 1:length(mid_x)){
      cd <- c(cd, convolute_and_distance(mid_x[[e]], mid_y[[e]], mid_z[[z]][[e]], f))
    }
    sums_cd <- c(sums_cd, sum(cd, na.rm = T))
  }
  
  return(list(g_select(sums_cd), z_index[which(sums_cd == g_select(sums_cd))[1]]))
}

#' @export
# convolutes two small MIDs and computes the distance (by f) between the convolution and the largest MID
convolute_and_distance <- function(mid_x, mid_y, mid_z, f){
  if (length(mid_x) < length(mid_y)) 
    return(f(convolute(mid_z, mid_x), mid_y)) else
      return(f(convolute(mid_z, mid_y), mid_x))
}

#' @export
# compute distance for each unique pair
calc_pair_distance <- function(pair, midata, f, g_select){
  print(pair)
  
  n_atoms_x <- midata$peak_n_atoms[pair[1]]
  mid_x <- lapply(1:length(midata$experiments), function(e) get_avg_mid(midata, pair[1], e))
  n_atoms_y <- midata$peak_n_atoms[pair[2]]
  mid_y <- lapply(1:length(midata$experiments), function(e) get_avg_mid(midata, pair[2], e))
  
  if (n_atoms_x == n_atoms_y){
    # compute distances for all experiments
    distance <- calc_distance_equal_c(mid_x, mid_y, f)
    result <- data.frame(metabolite_1 = pair[1], metabolite_2 = pair[2], 
                         middle_metabolite = distance[[2]], distance = distance[[1]])
    return(result)
  } 
  
  else {
    
    # number of carbons for the middle metabolite
    n_atoms_z <- abs(n_atoms_x - n_atoms_y)
    # indices of metabolites with n_atoms_z carbons
    z_index <- get_peak_index_n_atoms(midata, n_atoms_z)
    
    # make sure at least one metabolite with n_atoms_z carbons exists in the data
    if (is.null(z_index) == F) {
      # bring MIDs of these metabolites across all experiments
      mid_z <- lapply(z_index, function(p) lapply(1:length(midata$experiments), function(e) get_avg_mid(midata, p, e)))
      # compute the actual distances for each z, and select the one 
      distance <- calc_distance_unequal_c(mid_x, mid_y, mid_z, z_index, f, g_select)
      result <- data.frame(metabolite_1 = pair[1], metabolite_2 = pair[2], 
                           middle_metabolite = distance[[2]], distance = distance[[1]])
      return(result)
    } 
    # if not just return NA for this pair
    else return(data.frame(metabolite_1 = pair[1], metabolite_2 = pair[2], 
                           middle_metabolite = NA, distance = NA))
  }
  
}

#' @export
pairwise_matrix_v2 <- function(midata, f, g_select){
  
  # all unique pairs
  pairs <- combn(length(midata$peak_ids), 2)
  
  # compute distances and keep track of the convolutions
  results <- apply(pairs, 2, calc_pair_distance, midata, f, g_select)
  
  # convert dataframes into proper distance matrices
  distances <- sapply(results, function(x) return(x[4]))
  stopifnot(length(distances) == dim(pairs)[2])
  middle_mets <- sapply(results, function(x) return(x[3]))
  stopifnot(length(middle_mets) == dim(pairs)[2])
  
  dm <- matrix(NA, length(midata$peak_ids), length(midata$peak_ids))
  mmm <- matrix(NA, length(midata$peak_ids), length(midata$peak_ids))
  colnames(dm) <- rownames(dm) <- colnames(mmm) <- rownames(mmm) <- midata$peak_ids
  for (i in 1:ncol(pairs)){
    dm[pairs[1,i], pairs[2,i]] <- as.numeric(distances[i])
    dm[pairs[2,i], pairs[1,i]] <- as.numeric(distances[i])
    mmm[pairs[1,i], pairs[2,i]] <- as.numeric(middle_mets[i])
    mmm[pairs[2,i], pairs[1,i]] <- as.numeric(middle_mets[i])
  }
  
  result <- list(dm, mmm)
  names(result) <- c("distance_matrix", "middle_metabolite_matrix")
  
  return(result)
}

#' @export
# Define a function to validate if a matrix is a distance matrix
is_distance_matrix <- function(mat) {
  # necessary intervention to prevent the floating number error
  mat <- as.matrix(nearPD(mat)$mat)
  # Initialize a character vector to store failure reasons
  failures <- character(0)
  
  # Check if the input is a matrix
  if (!is.matrix(mat)) {
    failures <- append(failures, "Input is not a matrix")
  }
  
  # Check if the matrix is square
  if (nrow(mat) != ncol(mat)) {
    failures <- append(failures, "Matrix is not square")
  }
  
  # Check if the matrix is symmetric
  if (!isSymmetric(mat)) {
    failures <- append(failures, "Matrix is not symmetric")
  }
  
  # Check if the diagonal entries are zero
  if (!all(diag(mat) == 0)) {
    failures <- append(failures, "Diagonal entries are not zero")
  }
  
  # Check if the matrix is non-negative
  if (!all(mat >= 0)) {
    failures <- append(failures, "Matrix contains negative values")
  }
  
  # Check the triangle inequality
  for (i in 1:nrow(mat)) {
    for (j in 1:nrow(mat)) {
      for (k in 1:nrow(mat)) {
        if (round(mat[i,j], digits = 3) + round(mat[j,k], digits = 3) < round(mat[i,k], digits = 3)) {
          failures <- append(failures, "Triangle inequality fails")
        }
      }
    }
  }
  
  # If no failures, return NULL, else return the failure reasons
  if (length(failures) == 0) {
    return(NULL)
  } else {
    return(failures)
  }
}