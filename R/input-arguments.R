#
# functions for handling input parameters
#

#' Creates an InputData object to be used as input to similarity_matrix() or distance_matrix()
#' 
#' It takes as input a data frame with the required fields (see example below)
#' input_args <- data.frame(peak_areas_file_path = file.path(base_path, "Data", "simulation", "sim-data.tsv"),
#'    c13_correction = "F",
#'    measure = "cosine",
#'    fun = "cosine_sim",
#'    type = "similarity",
#'    perfection = "max",
#'    what_to_assign_to_na = 0,
#'    reaction_restriction = "F",
#'    reactions_file_path = "NULL",
#'    tolerance_ppm = 10,
#'    matrix_file_dir = paste0(local_file_path, type, "-matrices/"))
#'
#' @param input_args a data frame with fields listed above 
#' @export

parse_input_args <- function(input_args)
{
  # create an InputData structure
  input_data <- list()
  class(input_data) <- "InputData"
  
  # read peak areas from file if the input is a file path (character)
  # if (is.character(peak_areas)) 
  peak_areas <- as.data.frame(read.delim(input_args$peak_areas_file_path, header = T, sep = "\t"))
  
  if ("MassIsotopomer" %in% colnames(peak_areas))
    peak_areas <- peak_areas[, -which(colnames(peak_areas) == "MassIsotopomer")]
  
  # make an MIData object from peak_areas, fetching experiments from column names of the matrix
  midata <- MIData(peak_areas, colnames(peak_areas)[-(1:3)])
  
  # correct the MIDs if c13_correction is TRUE
  if (eval(parse(text = input_args$c13_correction)) == T){
    midata <- midata_transform(midata, c13correct)
    input_data$c13_correction <- "T"
  } else input_data$c13_correction <- "F"
    
  #
  input_data$midata <- midata
  #
  
  # similarity or distance function to apply to MIDs
  input_data$measure <- input_args$measure
  input_data$fun <- eval(parse(text = input_args$fun))
  input_data$type <- input_args$type
  input_data$perfection <- eval(parse(text = input_args$perfection))
  input_data$g_select <- eval(parse(text = input_args$g_select))
  input_data$get_middle_met_matrix <- eval(parse(text = input_args$get_middle_met_matrix))
  input_data$what_to_assign_to_na <- eval(parse(text = input_args$what_to_assign_to_na))
  #
  
  # add reaction restriction to input data
  if (input_args$reaction_restriction == "F" | input_args$reaction_restriction == F)
    input_data$reaction_restriction <- eval(parse(text = input_args$reaction_restriction)) else {
      input_data$reaction_restriction <- input_args$reaction_restriction
      input_data$reaction_data <- as.vector(read.delim(input_args$reactions_file_path, header = FALSE, sep = "\t")) 
    }
  
  
  # tolerance
  input_data$tolerance <- eval(parse(text = input_args$tolerance_ppm))*10^-6
  
  return(input_data)
}


#' @export
fetch_input_args <- function(input_row_index, input_file_name)
{
  options(warn = -1)
  input_args <- as.data.frame(read.delim(input_file_name, header = T, sep = "\t", check.names = F))
  return(input_args[input_row_index, ])
}
