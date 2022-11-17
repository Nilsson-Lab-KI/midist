#
# functions for handling input parameters
#

#' Creates an InputData object to be used as input to similarity_matrix()
#' 
#' It takes as input a data frame with the required fields (see example below)
#' input_args <- data.frame(peak_areas_file_path = file.path(base_path, "Data", "simulation", "sim-data.tsv"),
#'    c13_correction = "F",
#'    similarity_measure = "cosine_sim",
#'    reaction_restriction = "F",
#'    reactions_file_path = "NULL",
#'    tolerance_ppm = 10,
#'    similarity_matrix_file_dir = paste0(local_file_path, "similarity-matrices/"))
#'
#' @param input_args a dataframe with fields listed above 
#' @export

parse_input_args <- function(input_args)
{
  
  input_data <- list()
  class(input_data) <- "InputData"
  
  
  # make an MIData object from peak_areas
  peak_areas <- as.data.frame(read.delim(input_args$peak_areas_file_path, header = T, sep = "\t"))
  midata <- MIData(peak_areas, colnames(peak_areas)[-(1:3)])
  
  # correct the MIDs if c13_correction is TRUE
  if (eval(parse(text = input_args$c13_correction)) == T)
    midata <- midata_transform(midata, c13correct)
  
  #
  input_data$midata <- midata
  #
  
  # similarity function to apply to MIDs
  similarity <- similarity_measure <- eval(parse(text = input_args$similarity_measure))
  
  #
  input_data$similarity_name <- input_args$similarity_measure
  input_data$similarity <- similarity
  #
  
  # add reaction restriction to input data
  if (input_args$reaction_restriction == "F")
    input_data$reaction_restriction <- eval(parse(text = input_args$reaction_restriction)) else {
      input_data$reaction_restriction <- input_args$reaction_restriction
      input_data$reaction_data <- as.vector(read.delim(input_args$reactions_file_path, header = FALSE, sep = "\t")) 
    }
  
  
  # tolerance
  input_data$tolerance <- eval(parse(text = input_args$tolerance_ppm))*10^-6
  
  # c13 dir
  if (eval(parse(text = input_args$c13_correction)) == T)
    c13_dir <- paste0(input_args$similarity_matrix_file_dir, "with_13c_correction/") else
      c13_dir <- paste0(input_args$similarity_matrix_file_dir, "without_13c_correction/")
  
  # restriction sub-dir
  if (input_data$reaction_restriction == F)
    input_data$file_dir <- paste0(c13_dir, "no_restriction") else if (input_args$reaction_restriction == "formula")
      input_data$file_dir <- paste0(c13_dir, "formula_restriction") else
        input_data$file_dir <- paste0(c13_dir, "mass_restriction")
  
  return(input_data)
}



