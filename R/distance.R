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
apply_no_m0 <- function(f, x, y)
{
  return(f(x[-1], y[-1]))
}

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
  return(1 - cosine_sim(x,y))
}

#' @export
#'
cosine_dist_no_m0 <- function(x, y) apply_no_m0(cosine_dist, x, y)


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
  return(sqrt(sum(x*x)*sum(y*y)) - sum(x*y))
}

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
  diff <- x - y
  return(sum(diff*diff))
}

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
euclidean_dist <- function(x, y)
{
  return(sqrt(euclidean_dist_sq(x, y)))
}

#' @export
#'
euclidean_dist_no_m0 <- function(x, y) apply_no_m0(euclidean_dist, x, y)

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
    for(j in i:length(n_atoms)) {
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


#' Calculate similarity matrix, based on parameters specified in an InputData object
#' 
#' Calculates a similarity matrix based on the given similarity function for all peak pairs
#' for experiment e from an MIData object (midata).
#' If a list of reactions is input to input_data$reaction_data, the given similarity function will only be applied to peak pairs  
#' whose formula or mass difference (depending on "input_data$reaction_restriction") can be matched an allowed reaction
#'
#' @param e experiment index that corresponds to the same index in MIData
#' @param input_data InputData object (for more detail see ...)
#' @param write_to_file TRUE writes the similarity matrix to file. FALSE by default
#' @param return whether to return the similarity matrix. Choose FALSE for commandline runs
#' @export

pairwise_matrix_ds <- function(e, input_data, write_to_file = F, return = T)
{
  
  experiment <- input_data$midata$experiments[e]
  
  # print progress
  cat("Computing ", input_data$measure, " matrix for experiment: ",  experiment, '\n')
  if (input_data$reaction_restriction == F)
    cat("for all possible convolutions in data", '\n') else if (input_data$reaction_restriction == "mass")
      cat("for a restriction on mass difference", '\n') else
        cat("for a restriction on formula difference", '\n')
  
  
  # data dimensions
  n_metabolites <- length(input_data$midata$peak_ids)
  # met names
  met_names <- input_data$midata$peak_ids
  
  # create an empty matrix to be filled in with similarities
  pairwise_matrix <- matrix(NA, n_metabolites, n_metabolites)
  colnames(pairwise_matrix) <- rownames(pairwise_matrix) <- input_data$midata$peak_ids
  
  # loop over matrix elements (x,y) such that x <= y
  for (x in 1:n_metabolites) 
  {
    for (y in x:n_metabolites) 
    {
      
      if (input_data$reaction_restriction == F)
        pairwise_matrix[x,y] <- pairwise_matrix[y,x] <- conv_reduce(input_data$midata, 
                                                                    x, y, e, 
                                                                    input_data$fun, 
                                                                    input_data$perfection) 
      else # handle restrictions
      {
        reactions <- check_reactions(input_data)
        if (length(reactions) > 0)
          pairwise_matrix[x,y] <- pairwise_matrix[y,x] <- conv_reduce(input_data$midata, 
                                                                      x, y, e, 
                                                                      input_data$measure, 
                                                                      input_data$perfection, 
                                                                      input_data$what_to_assign_to_na)
      }
    }
  }
  
  
  if (write_to_file == T)
  {
    # # check if the directory exists and if not create it
    # print(paste0("Checking if ", input_data$file_dir, " exists, and if not creating it"))
    dir.create(input_data$file_dir, recursive = T)
    
    # write pairwise_matrix to file
    file_name <- file.path(input_data$file_dir, paste0(experiment, "_", input_data$measure, ".tsv"))
    write.table(pairwise_matrix, file_name, col.names = T, row.names = F, quote = F, sep = "\t")
  }
  
  if (return == T)
    return(pairwise_matrix)
  
}


#' @export
pairwise_matrix_all_ds <- function(e, input_data, write_to_file = F, return = T)
{
  
  experiment <- input_data$midata$experiments[e]
  
  # print progress
  cat("Computing ", input_data$measure, " matrix for experiment: ",  experiment, '\n')
  if (input_data$reaction_restriction == F)
    cat("for all possible convolutions in data", '\n') else if (input_data$reaction_restriction == "mass")
      cat("for a restriction on mass difference", '\n') else
        cat("for a restriction on formula difference", '\n')
  
  
  # data dimensions
  n_metabolites <- length(input_data$midata$peak_ids)
  # met names
  met_names <- input_data$midata$peak_ids
  
  # compute the pairwise matrix
  pairwise_matrix <- conv_reduce_all(input_data$midata, e, 
                                     input_data$fun, 
                                     input_data$perfection)
  colnames(pairwise_matrix) <- rownames(pairwise_matrix) <- input_data$midata$peak_ids
                  
  
  # # loop over matrix elements (x,y) such that x <= y
  # for (x in 1:n_metabolites) 
  # {
  #   for (y in x:n_metabolites) 
  #   {
  #     
  #     if (input_data$reaction_restriction == F)
  #       pairwise_matrix[x,y] <- pairwise_matrix[y,x] <- conv_reduce_all(input_data$midata, 
  #                                                                   e, 
  #                                                                   input_data$fun, 
  #                                                                   input_data$perfection) 
  #     else # handle restrictions
  #     {
  #       reactions <- check_reactions(input_data)
  #       if (length(reactions) > 0)
  #         pairwise_matrix[x,y] <- pairwise_matrix[y,x] <- conv_reduce_all(input_data$midata, 
  #                                                                     e, 
  #                                                                     input_data$measure, 
  #                                                                     input_data$perfection)
  #     }
  #   }
  # }
  
  
  if (write_to_file == T)
  {
    # # check if the directory exists and if not create it
    # print(paste0("Checking if ", input_data$file_dir, " exists, and if not creating it"))
    dir.create(input_data$file_dir, recursive = T)
    
    # write pairwise_matrix to file
    file_name <- file.path(input_data$file_dir, paste0(experiment, "_", input_data$measure, ".tsv"))
    write.table(pairwise_matrix, file_name, col.names = T, row.names = F, quote = F, sep = "\t")
  }
  
  if (return == T)
    return(pairwise_matrix)
  
}



check_reactions <- function(input_data, x, y)
{
  fun_main <- eval(parse(text = paste0("get_", input_data$reaction_restriction, "_difference")))
  fun <- eval(parse(text = paste0("get_", input_data$reaction_restriction)))
  
  # compute difference from the function above
  diff <- fun_main(list(fun(input_data$midata, x), fun(input_data$midata, y)))
  
  # compare this difference to the "allowed" reaction list
  # for formula
  if (is.character(diff))
  {
    reaction <- input_data$reaction_data[match(diff, input_data$reaction_data)]
    return(reaction)
  }
    
  
  
}


#' @export
combine_sqrt_sum <- function(pairwise_matrices){
  return(sqrt(Reduce('+', pairwise_matrices)))
}

#' @export
combine_sum <- function(pairwise_matrices){
  return(Reduce('+', pairwise_matrices))
}

#' @export
apply_summary_function <- function(vector, summary_function, exclude_missing_values = T)
{
  
  new_vector <- vector[which(is.na(vector) == F & is.infinite(vector) == F)]
  
  if (exclude_missing_values == F)
    new_vector <- vector
  
  return(summary_function(new_vector))
}

#' @export
combine_max <- function(pairwise_matrices)
{
  return(matrix(apply(do.call(rbind.data.frame, 
                              lapply(pairwise_matrices, as.vector)), 
                      2, 
                      apply_summary_function, 
                      max), 
                ncol = ncol(pairwise_matrices[[1]]), 
                byrow = F))
}

#' @export
combine_min <- function(pairwise_matrices)
{
  return(matrix(apply(do.call(rbind.data.frame, 
                              lapply(pairwise_matrices, as.vector)), 
                      2, 
                      apply_summary_function, 
                      min), 
                ncol = ncol(pairwise_matrices[[1]]), 
                byrow = F))
}

#' @export
combine_mean <- function(pairwise_matrices)
{
  return(matrix(apply(do.call(rbind.data.frame, 
                              lapply(pairwise_matrices, as.vector)), 
                      2, 
                      apply_summary_function, 
                      mean), 
                ncol = ncol(pairwise_matrices[[1]]), 
                byrow = F))
}

#' @export
combine_median <- function(pairwise_matrices)
{
  return(matrix(apply(do.call(rbind.data.frame, 
                              lapply(pairwise_matrices, as.vector)), 
                      2, 
                      apply_summary_function, 
                      median), 
                ncol = ncol(pairwise_matrices[[1]]), 
                byrow = F))
}

