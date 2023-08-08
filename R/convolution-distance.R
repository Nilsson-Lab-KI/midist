#
# Convolution-based MID distance measures
#


#' Convolution-based measure between MIDs for peak x and y for
#' experiment e from an MIData object (mi_data).
#'
#' If x, y have the same number of atoms, f(x,y) is returned.
#' If x is smaller than y, this function computes f(x*z, y) for all
#' peaks z such that x*z is the same size as y (and vice versa), and calls g()
#' on the resulting vector of values from f to select a value.
#' For example, setting g = which.max gives the maximum value.
#'
#' @param mi_data the MIdata object
#' @param x peak index
#' @param y peak index
#' @param e experiment index
#' @param f A function f(x, y) taking two MIDs and returning a scalar.
#' @param g a function g taking a non-empty vector of values and returning
#' the index of the "best" element; for example, which.min
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
  # find the MIDs of metabolite with index x, y
  mid_x <- get_avg_mid(mi_data, x, e)
  mid_y <- get_avg_mid(mi_data, y, e)

  if (n_atom_x == n_atom_y) {
    # equal numbers just calculate f
    return(list(values = f(mid_x, mid_y), index = NA))
  }
  else {
    # x is strictly smaller than y
    # get MIDs of metabolites z to convolute with
    n_atom_z <- n_atom_y - n_atom_x
    z_index <- get_peak_index_n_atoms(mi_data, n_atom_z)
    if (length(z_index) > 0) {
      # this is either an MI x z_index matrix if e is scalar,
      # or an MI x experiments x z_index array if e is a vector
      mids_z <- sapply(z_index,
                       function(i) get_avg_mid(mi_data, i, e),
                       simplify = "array")
      # compute all convolutions x*z for each z
      if(length(e) == 1)
        mids_xz <- convolute_cols(mid_x, mids_z)
      else {
        #cat("mid_x dim = ", dim(mid_x))
        #cat("mids_z dim = ", dim(mids_z))
        # this yields an MI x experiments x z_index array
        mids_xz <- convolute_array(mid_x, mids_z)
      }
      # calculate f between y and all x*y and add indices
      f_values <- apply(mids_xz, MARGIN = 2, f, mid_y)
      # return best value and index
      f_index <- g(f_values)
      return(list(values = f_values[f_index], index = z_index[f_index]))
    } else {
      # no matching metabolites to convolute with
      return(list(values = NA, index = NA))
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
#' @param e Experiment index. If e is a vector, the function f(x, y)
#' must accept matrices.
#' @param f A function f(x, y) taking two MIDs, or two matrices whose columns are MIDs.
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
    n_atoms_x <- n_atoms[i]
    # get all MIDs for x
    x_index <- get_peak_index_n_atoms(mi_data, n_atoms_x)

    # for each metabolites size for y, larger than x
    for(j in i:length(n_atoms)) {
      n_atoms_y <- n_atoms[j]
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
      # copy block to full matrices
      conv_values[x_index, y_index] <- block$values
      conv_values[y_index, x_index] <- t(block$values)
      conv_index[x_index, y_index] <- block$index
      conv_index[y_index, x_index] <- t(block$index)
    }
  }
  return(list(value = conv_values, index = conv_index))
}


# apply a selection function g and return a list of
# g values and corresponding indices
g_list <- function(mids_y, mids_xz, z_index, f, g)
{
  n_y <- dim(mids_y)[2]
  values <- rep(as.double(NA), n_y)
  index <- rep(as.integer(NA), n_y)
  for(i in 1:n_y) {
    f_values <- apply(mids_xz, MARGIN = 2, f, mids_y[, i])
    f_index <- g(f_values)
    # store values and index separately
    values[i] <- f_values[f_index]
    index[i] <- z_index[f_index]
  }
  # add index vector as attribute
  return(list(values = values, index = index))
}



#' compute conv_reduce for one matrix block where all x are the same size
#' and all y are the same size, but y is larger than x
#'
#' @param mi_data An MIData object
#' @param e experiment index
#' @param f A distance function
#' @param g A selection function, such as min_nonempty
#' @param x_index index of the peak x
#' @param y_index index of the peak y
#' @param z_index index of the peak z
#'
conv_reduce_block <- function(mi_data, e, f, g, x_index, y_index, z_index)
{
  n_x <- length(x_index)
  n_y <- length(y_index)
  ## allocate matrices
  block = list(values = matrix(as.double(NA), n_x, n_y),
               index = matrix(as.integer(NA), n_x, n_y))
  if(length(z_index) > 0) {
    # get MID matrices, each column an MID
    mids_x <- sapply(x_index, function(i) get_avg_mid(mi_data, i, e))
    mids_y <- sapply(y_index, function(i) get_avg_mid(mi_data, i, e))
    mids_z <- sapply(z_index, function(i) get_avg_mid(mi_data, i, e))

    for(i in 1:n_x) {
      # compute all convolutions x*z for each z
      mids_xz <- convolute_cols(mids_x[, i], mids_z)
      # calculate g(f(x*z, y) ...) for all y (rows) and all convolutions x*z
      g_result <- g_list(mids_y, mids_xz, z_index, f, g)
      # NOTE: is.vector(x) is FALSE when x is a vector but has attributes!
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
  mids_y <- sapply(y_index, function(i) get_avg_mid(mi_data, i, e))

  # calculate f(x,y) for each x,y (no attributes in this case)
  for(i in 1:n_x) {
    block$values[i, ] <- apply(mids_y, MARGIN = 2,
                               function(mid_y) f(mids_x[, i], mid_y))
  }
  return(block)
}


#' Filter a matrix so that a given fraction of its real-valued elements
#' are set to the corresponding quantile
#' @param pairwise_matrix A symmetric matrix
#' @param percentile The fraction (NOT a percentile) of elements to keep
filter_pairwise_matrix_global <- function(pairwise_matrix, percentile = 0.01)
{
  # remove diagonals
  diag(pairwise_matrix) <- NA

  # create an empty matrix to be filled in only by those who pass the filtering criteria
  filtered_pm <- matrix(NA, nrow(pairwise_matrix), ncol(pairwise_matrix))
  colnames(filtered_pm) <- rownames(filtered_pm) <- colnames(pairwise_matrix)

  # get indices of non-infinite, non-NA elements
  non_inf_ind <- which(is.infinite(pairwise_matrix) == F & is.na(pairwise_matrix) == F)
  # exclude infinites and NAs
  non_inf_vec <- pairwise_matrix[non_inf_ind]
  # compute the percentile (top 1%)
  threshold <- as.numeric(stats::quantile(non_inf_vec, probs = percentile))
  # fill in
  filtered_pm[as.numeric(non_inf_ind[which(non_inf_vec <= threshold)])] <-
    pairwise_matrix[as.numeric(non_inf_ind[which(non_inf_vec <= threshold)])]
  return(filtered_pm)
}


#' Returns the metabolite index which was chosen by fun.
#' It is supposed to be used in the combine() function.
#' @param vector a row from a pairwise matrix
#' @param fun a function selecting an element of the vector, like max() or min()
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
#' @param fun Function to apply
#' @returns a list of 3 matrices: (1) the combined pairwise matrix, (2) the combined
#' middle metabolite matrix, and (3) the experiment index chosen for each element
#' @export
combine <- function(pairwise_matrices, middle_met_matrices, fun)
{
  # flatten matrices
  pairwise_vec <- do.call(rbind.data.frame, lapply(pairwise_matrices, as.vector))
  middle_met_vec <- do.call(rbind.data.frame, lapply(middle_met_matrices, as.vector))
  # find the experiment
  index <- as.vector(unlist(apply(pairwise_vec, 2, get_fun_index, fun)))
  # combined pairwise matrix
  pm <- matrix(
    unlist(lapply(1:length(index),
                  function(x, pairwise_vec, index) pairwise_vec[[x]][index[[x]][1]],
                  pairwise_vec, index)),
    ncol = ncol(pairwise_matrices[[1]]),
    byrow = F
  )
  # combined middle metabolite matrix
  mmm <- matrix(
    unlist(lapply(1:length(index),
                  function(x, pairwise_vec, index) middle_met_vec[[x]][index[[x]][1]],
                  vector, index)),
    ncol = ncol(pairwise_matrices[[1]]),
    byrow = F
  )
  # matrix of experiment indices chosen by fun
  ei <- matrix(index,
    ncol = ncol(pairwise_matrices[[1]]),
    byrow = F
  )

  return(list(pm, mmm, ei))
}


# compute distances for equal carbon metabolites, and make a selection based on g_select
calc_distance_equal_c <- function(mid_x, mid_y, f){

  cd <- c()
  for (e in 1:length(mid_x)){
    cd <- c(cd, f(mid_x[[e]], mid_y[[e]]))
  }
  # returns 2 elements: (1) the final distance, (2) the middle metabolite which was chosen by g_select
  return(list(sum(cd, na.rm = T), NA))
}


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


# convolutes two small MIDs and computes the distance (by f) between the convolution and the largest MID
convolute_and_distance <- function(mid_x, mid_y, mid_z, f){
  if (length(mid_x) < length(mid_y))
    return(f(convolute(mid_z, mid_x), mid_y)) else
      return(f(convolute(mid_z, mid_y), mid_x))
}


# compute distance for each unique pair
calc_pair_distance <- function(pair, midata, f, g_select)
{
#  print(pair)
  n_atoms_x <- midata$peak_n_atoms[pair[1]]
  mid_x <- lapply(1:length(midata$experiments), function(e) get_avg_mid(midata, pair[1], e))
  n_atoms_y <- midata$peak_n_atoms[pair[2]]
  mid_y <- lapply(1:length(midata$experiments), function(e) get_avg_mid(midata, pair[2], e))

  if (n_atoms_x == n_atoms_y) {
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


# Compute distance and middle metabolite matrices
# Returns a list of length 2
pairwise_matrix_v2 <- function(midata, f, g_select)
{
  # all unique pairs
  pairs <- utils::combn(length(midata$peak_ids), 2)

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


#' Compute a pairwise distance matrix and optionally save it
#' @param midata An MIData object
#' @param f A distance function
#' @param g_select A selection function
#' @param rdata_fname A file name to save results in, if return == F
#' @param return A boolean; if FALSE, results are saved using base::save
#' @returns A pairwise matrix if return == TRUE; else nothing
#' @export
remn_v2 <- function(midata, f, g_select, rdata_fname, return = T)
{
  remn_output <- pairwise_matrix_v2(midata, f, g_select)

  if (return == T)
    return(remn_output)
  else {
    distance_matrix <- remn_output$distance_matrix
    middle_metabolite_matrix <- remn_output$middle_metabolite_matrix
    save(midata,
         distance_matrix,
         middle_metabolite_matrix,
         f,
         g_select,
         file = rdata_fname)
  }
}


#' Compute distance matrix from isotopic enrichment
#' @param midata An MIData object
#' @param experiments A list of experiments (names) to use
#' @param method Name of  distance metric, as in stats::dist. Use either
#' of "euclidean", "manhattan", or "canberra"
#' @export
enrichment_dist_matrix <- function(midata, experiments, method = "euclidean")
{
  isotopic_enrichments <- do.call(
    rbind.data.frame,
    lapply(1:length(midata$peak_ids),
           function(p) apply(as.matrix(get_avg_mid(midata, p)), 2, isotopic_enrichment)
    )
  )
  # compute distance matrix
  enrichments <- as.matrix(
    stats::dist(
      isotopic_enrichments[, midata$experiments %in% experiments],
      method = method
    )
  )
  rownames(enrichments) <- colnames(enrichments) <- midata$peak_ids
  return(enrichments)
}


#' Compute pairwise distance matrices for individual experiments
#'
#' @param midata An MIData object
#' @param f A distance function, as in conv_reduce
#' @param g_select A function for selecting best convolutions, as in conv_reduce
#' @param rdata_fname NOT USED
#' @param return NOT USED
#' @returns A list with elements $value holding a list of distance matrices,
#' and $index holding a list of "middle metabolite" index matrices
#' @export
remn_v1 <- function(midata, f, g_select, rdata_fname, return = T)
{
  values <- list()
  index <- list()
  for (e in 1:length(midata$experiments)) {
    assign_list[values[[e]], index[[e]]] <- conv_reduce_all(midata, e, f, g_select)
  }
  return(list(values = values, index = index))
}

