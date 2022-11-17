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

calc_dot <- function(x, y) {1 - sum(x*y)}

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

cosine_sim_no_m0 <- function(x, y) apply_no_m0(similarity_measure, x, y)


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

similarity_matrix <- function(e, input_data, write_to_file = F, return = T)
{
  
  experiment <- input_data$midata$experiments[e]
  cat("Computing ", input_data$similarity_name, " matrix for experiment ",  experiment, '\n')
  if (input_data$reaction_restriction == F)
    cat("for all possible convolutions in data") else if (input_data$reaction_restriction == "mass")
      cat("for a restriction on mass difference") else
        cat("for a restriction on formula difference")
  
  
  # data dimensions
  n_metabolites <- length(input_data$midata$peak_ids)
  
  # met names
  met_names <- input_data$midata$peak_ids
  
  # create an empty matrix to be filled in with similarities
  sm <- matrix(NA, n_metabolites, n_metabolites)
  colnames(sm) <- input_data$midata$peak_ids
  rownames(sm) <- input_data$midata$peak_ids
  
  # loop over matrix elements (x,y) such that x <= y
  for (x in 1:n_metabolites) 
  {
    for (y in x:n_metabolites) 
    {
      # if type == "mass"
      if (input_data$reaction_restriction == "mass")
      {
        # check if the mass difference of the two metabolites match an allowed reaction
        mass_diff <- abs(get_mass(input_data$midata, x) - get_mass(input_data$midata, y))
        if (length(compare_mass_diff_to_list(mass_diff, input_data$reaction_data, input_data$tolerance)) > 0)
          sm[x,y] <- sm[y,x] <- conv_similarity(input_data$midata, x, y, e, input_data$similarity)
      }
      # if type == "formula"
      else if (input_data$reaction_restriction == "formula")
      {
        # check if the formula difference of the two metabolites match an allowed reaction
        formula_diff <- get_formula_difference(list(get_formula(input_data$midata, x), get_formula(input_data$midata, y)))
        if (length(which(is.na(unlist(lapply(formula_diff, match, input_data$reaction_data))) != T)) > 0)
          sm[x,y] <- sm[y,x] <- conv_similarity(input_data$midata, x, y, e, input_data$similarity)
      }
      else
      {
        sm[x,y] <- sm[y,x] <- conv_similarity(input_data$midata, x, y, e, input_data$similarity)
      }
      
    }
  }
  
  if (write_to_file == T)
  {
    # check if the directory exists and if not create it
    print(paste0("Checking if ", input_data$file_dir, " exists, and if not creating it"))
    dir.create(input_data$file_dir)
    
    # write sm to file
    file_name <- file.path(input_data$file_dir, paste0(experiment, "_", input_data$similarity_name, ".tsv"))
    write.table(sm, file_name, col.names = T, row.names = F, quote = F, sep = "\t")
    cat("'", input_data$similarity_name, "'", " matrix for experiment '", experiment, "' was written to file:", '\n',  file_name)
  }
  
  if (return == T)
    return(sm)
  
}

