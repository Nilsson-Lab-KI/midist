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
#' @param similarity_measure A similarity function on MIDs
#' @param remove_m0 whether to remove M+0 from MIDs (TRUE if MIDs are corrected for the natural 13C. FALSE by default)
#' @returns the resulting similarity matrix
#' @export

similarity_matrix <- function(midata, e, similarity_measure, remove_m0 = FALSE)
{
  
  # c13 correct midata, and update similarity to remove M+0
  if (remove_m0 == TRUE)
  {
    midata <- midata_transform(midata, c13correct)
    similarity <- function(x,y) apply_no_m0(similarity_measure, x,y)
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


#' Calculate similarity matrix, IF THE FORMULA OR MASS DIFFERENCE IS ACCEPTABLE, 
#' based on a specific similarity/distance measure, for metabolites from an MIData object
#'
#' Calculates a similarity matrix based on the given similarity function for all peak pairs  
#' whose formula or mass difference can be matched an allowed reaction
#' for experiment e from an MIData object (midata).
#'
#' @param midata the MIdata object
#' @param e experiment index
#' @param similarity_measure A similarity function on MIDs
#' @param allowed_reactions a vector of allowed reactions. Either formulas or masses depending on whether type is "formulas" or "masses", respectively
#' @param type how to asses the intermediate reaction - based on whether formula or mass.
#' @param remove_m0 whether to remove M+0 from MIDs (TRUE if MIDs are corrected for the natural 13C. FALSE by default)
#' @returns the resulting similarity matrix
#' @export

similarity_matrix_with_reaction_filter <- function(midata, e, similarity_measure, allowed_reactions, type = "formula", 
                                                   remove_m0 = FALSE, tol = 10)
{
  
  # c13 correct midata, and update similarity to remove M+0
  if (remove_m0 == TRUE)
  {
    midata <- midata_transform(midata, c13correct)
    similarity <- function(x,y) apply_no_m0(similarity_measure, x,y)
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
      
      # check if the formula difference of the two metabolites match an allowed reaction
      if (type == "mass"){
        mass_diff <- abs(get_mass(midata, x) - get_mass(midata, y))
        if (length(compare_mass_diff_to_list(mass_diff, allowed_reactions, tol)) > 0)
          sm[x,y] <- sm[y,x] <- conv_similarity(midata, x, y, e, similarity)
      }
      # if type == "formula"
      else {
        formula_diff <- get_formula_difference(list(get_formula(midata, x), get_formula(midata, y)))
        if (length(which(is.na(unlist(lapply(formula_diff, match, allowed_reactions))) != T)) > 0)
          sm[x,y] <- sm[y,x] <- conv_similarity(midata, x, y, e, similarity)
      }
      
    }
  }
  
  return(sm)
  
}


########## DON'T REALLY KNOW WHERE TO PUT THESE FUNCTIONS ###########


# formulas is a list of length two, that stores the formulas to be compared to each other
# the output is a list that represents the reaction in both directions
#' @export
get_formula_difference <- function(formulas)
{
  broken_formulas <- lapply(formulas, break_formula)
  diff <- list()
  diff[[1]] <- data.frame(Symbol = broken_formulas[[1]]$Symbol,
                          ElementNumber = broken_formulas[[1]]$ElementNumber - broken_formulas[[2]]$ElementNumber)
  diff[[2]] <- data.frame(Symbol = broken_formulas[[1]]$Symbol,
                          ElementNumber = broken_formulas[[2]]$ElementNumber - broken_formulas[[1]]$ElementNumber)
  
  return(lapply(diff, merge_formula))
}
#' @export
break_formula <- function(formula)
{
  # get the periodic table elements (symbols)
  Symbol <- c("C", "H", "He", "Li", "Be", "B", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu",
              "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce",
              "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr",
              "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv",
              "Ts", "Og")
  
  # make a data frame for the input formula
  brokenFormula <- as.data.frame(Symbol)
  brokenFormula$ElementNumber <- 0
  
  # splitting the formula into the smallest pieces
  formula_split <- unlist(strsplit(formula, ""))
  
  # merging the lower case letters (if there is any) to their upper case letter
  grep_lower <- grepl("[a-z]", formula_split)
  l <- which(grep_lower == TRUE)
  if (length(l) > 0){
    formula_split[l-1] <- paste(formula_split[l-1], formula_split[l], sep = "")
    formula_split[l] <- NA
    formula_split <- formula_split[-which(is.na(formula_split)==TRUE)]
  }
  
  
  # find elements' location in the broken formula
  element_ind <- which(formula_split %in% Symbol)
  
  # now we need to assign "1" in case the last element doesn't have a number
  if (element_ind[length(element_ind)] == length(formula_split)){
    formula_split <- c(formula_split, "1")
  }
  
  # now we need to assign "1" to any element in the middle that doesn't have a number
  if (length(element_ind) > 1){
    for (i in 1:(length(element_ind)-1)){
      row_index <- which(brokenFormula$Symbol == formula_split[element_ind[i]])
      if (element_ind[i]+1 == element_ind[i+1]){
        brokenFormula$ElementNumber[row_index] <- 1
      }
      if (element_ind[i]+1 != element_ind[i+1]){
        brokenFormula$ElementNumber[row_index] <- as.numeric(paste(formula_split[(element_ind[i]+1):(element_ind[i+1]-1)], collapse = ""))
      }
    }
  }
  
  i <- length(element_ind)
  row_index <- which(brokenFormula$Symbol == formula_split[element_ind[i]])
  brokenFormula$ElementNumber[row_index] <- as.numeric(paste(formula_split[(element_ind[i]+1):length(formula_split)], collapse = ""))
  # brokenFormula <- brokenFormula[which(brokenFormula$ElementNumber!=0),]
  return(brokenFormula)
}
#' @export
get_formula_mass <- function(formula, periodic_table)
{
  broken_formula <- break_formula(formula)
  broken_formula$mass <- periodic_table$AtomicMass[match(broken_formula$Symbol, periodic_table$Symbol)] * broken_formula$ElementNumber
  return(sum(broken_formula$mass))
}
#' @export
get_mass <- function(midata, p){
  return(
    midata$peak_masses[p])
}
#' @export
merge_formula <- function(broken_formula)
{
  broken_formula <- broken_formula[which(broken_formula$ElementNumber != 0),]
  merged_formula <- as.vector(apply(as.data.frame(apply(broken_formula, 1, paste0, collapse = "")), 2, paste0, collapse = ""))
  return(gsub(" ", "", merged_formula))
}
#' @export
compare_formula_diff_to_list <- function(formula_diff, allowed_reactions)
{
  ind <- which(is.na(unlist(lapply(formula_diff, match, allowed_reactions))) != T)
  if (length(ind) > 0){
    return(formula_diff[[ind]])
  }
}
#' @export
compare_mass_diff_to_list <- function(mass_diff, allowed_masses, tol = 10)
{
  error_l <- mass_diff - mass_diff*tol*10^-6
  error_u <- mass_diff + mass_diff*tol*10^-6
  
  ind <- which(allowed_masses > error_l & allowed_masses < error_u)
  if (length(ind) > 0){
    return(allowed_masses[ind])
  }
}
#' @export
get_formula <- function(mi_data, p)
{
  return(
    mi_data$peak_formulas[p])
}
