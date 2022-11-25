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


#' @export
#'
cosine_sim_no_m0 <- function(x, y) apply_no_m0(cosine_sim, x, y)

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
  conv_reduce(mi_data, x, y, e, similarity, max, 0)
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
conv_distance <- function(mi_data, x, y, e, distance, default = Inf)
{
  conv_reduce(mi_data, x, y, e, distance, min, default)
}


#' Calculates f(x, y) between MIDs for peak x and y for
#' experiment e from an MIData object (mi_data)
#' If x,y have unequal atom numbers, computes f(x*z, y) for all atoms of the same
#' size as y (if y is the longer vector) and returns g( ... ) of the list of
#' value of f. For example f = max gives the maximum value.
#' if similarity returns NA for some pair, the NA value is ignored when taking
#' maximum similarity; if there are no possible convolutions, returns 0
#'
#' @param mi_data the MIdata object
#' @param x peak index
#' @param y peak index
#' @param e experiment index
#' @param f A function f(x, y) taking two MIDs.
#' @param g a function g taking a vector of values f1, f2, ...
#' @param default value to use when there is no possible convolution
#' @returns the resulting value g(f1, f2, ...)
#' @export

conv_reduce <- function(mi_data, x, y, e, f, g, default = 0)
{
  # number of carbon atoms of metabolites x, y
  n_atom_x <- get_peak_n_atoms(mi_data, x)
  n_atom_y <- get_peak_n_atoms(mi_data, y)
  # ensure MID x is smaller or equal to MID y
  if(n_atom_x > n_atom_y) {
    return(conv_reduce(mi_data, y, x, e, f, g, default))
  }

  # find the MIDs of metabolite with index x, y from experiment e
  mid_x <- get_avg_mid(mi_data, x, e)
  mid_y <- get_avg_mid(mi_data, y, e)

  if (n_atom_x == n_atom_y) {
    # in case of equal number, just calculate f
    return(f(mid_x, mid_y))
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

      # calculate f between the larger metabolite and all possible convolutions
      f_values <- unlist(lapply(convolutions, f, mid_y))
      # remove any missing values
      f_values <- f_values[!is.na(f_values)]
      if(length(f_values) == 0) {
        # all NA, so return NA
        return(NA)
      }
      else {
        # return the function g
        return(g(f_values))
      }
    }
    else {
      # no matching metabolites to convolute with
      # return default
      return(default)
    }
  }
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
                                                                    input_data$perfection, 
                                                                    input_data$what_to_assign_to_na) 
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
    # check if the directory exists and if not create it
    print(paste0("Checking if ", input_data$file_dir, " exists, and if not creating it"))
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
combine_max <- function(pairwise_matrices){
  return(Reduce(max, pairwise_matrices))
}

#' @export
combine_min <- function(pairwise_matrices){
  return(Reduce(min, pairwise_matrices))
}

#' @export
combine_mean <- function(pairwise_matrices){
  return(Reduce(mean, pairwise_matrices))
}

#' @export
combine_median <- function(pairwise_matrices){
  return(Reduce(median, pairwise_matrices))
}
