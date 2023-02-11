#
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
    return(NA)
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

#' Calculates a convolution-based measure between MIDs for peak x and y for
#' experiment e from an MIData object (mi_data).
#'
#' The function f(x, y) must return a real value or NA for any x, y with the
#' same number of atoms.
#' If x, y have the same number of atoms, f(x,y) is returned.
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

conv_reduce <- function(mi_data, x, y, e, f, g) {
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
  } else {
    # x is strictly smaller than y
    # get MIDs of metabolites z to convolute with
    n_atom_z <- n_atom_y - n_atom_x
    mids_z <- get_avg_mids_by_size(mi_data, n_atom_z, e)
    if (length(mids_z) > 0) {
      # compute all convolutions x*z for each z
      convolutions <- convolute_cols(mid_x, mids_z)
      # calculate f between the larger metabolite and all possible convolutions
      f_values <- apply(convolutions, MARGIN = 2, f, mid_y)
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
  # allocate square matrix
  n_met <- length(mi_data$peak_ids)
  result <- matrix(0, nrow = n_met, ncol = n_met)

  # unique metabolite sizes, sorted
  n_atoms <- sort(unique(mi_data$peak_n_atoms))
  # for each metabolite size for x
  for(i in 1:length(n_atoms)) {
    n_atoms_x <- n_atoms[[i]]
    x_index <- get_peak_index_n_atoms(mi_data, n_atoms_x)

    # for each metabolites size for y, larger than x
    for(j in i:length(n_atoms)) {
      n_atoms_y <- n_atoms[[j]]
      # get all MIDs y of this size
      y_index <- get_peak_index_n_atoms(mi_data, n_atoms_y)
      # matrix, each column one MID
      mids_y <- sapply(y_index, function(i) get_avg_mid(mi_data, i, e))
      stopifnot(is.matrix(mids_y))

      if(n_atoms_x == n_atoms_y) {
        # if y has same size as x, just calculate f(x,y)
        for(xi in x_index) {
          mid_x <- get_avg_mid(mi_data, xi, e)
          result[xi, y_index] <- result[y_index, xi] <-
            apply(mids_y, MARGIN = 2, function(y) g(c(f(mid_x, y))))
        }
      }
      else {
        # metabolites z to convolute with each x
        n_atoms_z <- n_atoms_y - n_atoms_x
        z_index <- get_peak_index_n_atoms(mi_data, n_atoms_z)
        if(length(z_index) > 0) {
          # matrix, each column one MID
          mids_z <- get_avg_mids_by_size(mi_data, n_atoms_z, e)
          if(!is.matrix(mids_z)) {
            print(z_index)
            print(mids_z)
            stop()
          }
          for(xi in x_index) {
            # MID of metabolite x
            mid_x <- get_avg_mid(mi_data, xi, e)
            # compute all convolutions x*z for each z
            xz_mids <- convolute_cols(mid_x, mids_z)
            stopifnot(is.matrix(xz_mids))
            # calculate f(x*z, y) for all y (rows) and all convolutions x*z
            # apply g() here?
            f_values <- apply(mids_y, MARGIN = 2,
                              function(y) g(apply(xz_mids, MARGIN = 2, f, y)))
            if(!is.vector(f_values)) {
              print(f_values)
              stop()
            }
            result[xi, y_index] <- result[y_index, xi] <- f_values

            #if(length(z_index) == 1){
            #  result[xi, y_index] <- result[y_index, xi] <- f_values
            #}
            #else {
            #  # compute g(f(y, ..), ...) for each y and store
            #  result[xi, y_index] <- result[y_index, xi] <- apply(f_values, 2, g)
            #}
          }
        }
        else {
          # no metabolites z to convolute x with
          for(xi in x_index) {
            result[xi, y_index] <- result[y_index, xi] <-
              rep(suppressWarnings(g(c())), n = length(y_index))
          }
        }
      }
    }
  }
  return(result)
}

#' Same as conv_reduce_all, but uses a function g_select to pick a particular
#' convolution, and also returns the matrix of peak indices picked by g_select
#' For example g = max gives the maximum value, and g_select = which.max gives
#' the index of the maximum value.
#' The function g must handle empty lists g(c()) in case there are no possible convolutions.
#' In that case, g_select will not be applied.
#'
#' @param mi_data the MIdata object
#' @param e experiment index
#' @param f A function f(x, y) taking two MIDs.
#' @param g a function g taking a vector of values f1, f2, ...
#' @param g_select either which.max() or which.min() dependent on g being max() or min()
#' @returns the matrix of g(f(x,y) ...) and g_select(f(x,y) ...) values for all x,y
#' @export
conv_reduce_all_select <- function(mi_data, e, f, g, g_select)
{
  # allocate square matrix
  n_met <- length(mi_data$peak_ids)
  result <- matrix(0, nrow = n_met, ncol = n_met)
  middle_met_matrix <- matrix(NA, nrow = n_met, ncol = n_met)

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
      mids_z <- get_avg_mids_by_size(mi_data, n_atoms_z, e)

      for(x in get_peak_index_n_atoms(mi_data, n_atoms_x)) {
        # MIDs of metabolite with index x
        mid_x <- get_avg_mid(mi_data, x, e)

        if(n_atoms_x == n_atoms_y) {
          # for metabolites y of same size as x, just calculate f(x,y)
          result[x, y_index] <- result[y_index, x] <- apply(mids_y,
            MARGIN = 2,
            function(y) g(c(f(mid_x, y)))
          )
        } else {
          if (length(mids_z) > 0) {
            # get indices of all z
            z_ind <- get_peak_index_n_atoms(mi_data, n_atoms_z)
            # compute all convolutions x*z for each z
            convolutions <- convolute_cols(mid_x, mids_z)
            # for each y, calculate g(f(x*z, y)) over all convolutions x*z
            # get pairwise measures for each convolution
            pm_conv <- apply(mids_y, 2, function(y) apply(convolutions, MARGIN = 2, f, y))


            if (length(z_ind) == 1){
              result[x, y_index] <- result[y_index, x] <- pm_conv
              middle_met_matrix[x, y_index] <- middle_met_matrix[y_index, x] <- z_ind
            }
            else {
              result[x, y_index] <- result[y_index, x] <- apply(pm_conv, 2, g)
              middle_met_matrix[x, y_index] <- middle_met_matrix[y_index, x] <- z_ind[apply(pm_conv, 2, g_select)]
            }

          }
          else {
            # no metabolites z to convolute x with
            result[x, y_index] <- result[y_index, x] <- rep(suppressWarnings(g(c())), n = length(y_index))
            middle_met_matrix[x, y_index] <- middle_met_matrix[y_index, x] <- NA

          }
        }
      }
    }
  }
 return(list(result, middle_met_matrix))
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
subsample_pairwise_matrix <- function(pairwise_matrix, percentile = 0.01){
  
  
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
combine_sqrt_sum <- function(pairwise_matrices) {
  return(sqrt(Reduce("+", pairwise_matrices)))
}

#' Combines pairwise matrices across multiple experiments by their sum.
#' Note that this cannot be used if the middle metabolite matrix and experiment index matrix are to be tracked.
#'
#' @param pairwise_matrices a list of pairwise matrices across different experiments
#' @export
combine_sum <- function(pairwise_matrices) {
  return(Reduce("+", pairwise_matrices))
}


#' Returns the metabolite index which was chosen by fun. It is supposed to be used in the combine() function.
#' @param vector a row from a pairwise matrix
#' @param fun a function like max() or min()
#' @export
get_fun_index <- function(vector, fun) {
  vector <- as.vector(vector)
  non_na_inf_ind <- which(is.na(vector) == F & is.infinite(vector) == F)
  new_vector <- vector[non_na_inf_ind]

  return(non_na_inf_ind[which(new_vector == fun(new_vector))][1])
}


#' Combines pairwise matrices across multiple experiments based on the fun function.
#' Suitable for cases where the middle metabolite matrix is known and experiment index matrix is to be tracked.
#'
#' @param pairwise_matrices a list of pairwise matrices across different experiments
#' @param middle_met_matrices a list of middle metabolite matrices across different experiments, where a middle metabolite denotes the metabolite that was chosen as the best convolution for a given metabolite pair of unequal carbon number
#'
#' @export
combine <- function(pairwise_matrices, middle_met_matrices, fun)
{
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
