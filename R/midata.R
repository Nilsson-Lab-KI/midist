#
# MID data functions
#

#' Construct an MIData (mass isotopomer data) object
#'
#' @param peak_areas a peak_area data.frame. The first column of peak_areas must be peak identifiers,
#' repeated for each MI of the same peak, and MIs must be increasing 0,1,...,n for each peak
#' @param exp_names a optional list of experiment names matching columns 2,3... in peak_areas
#' @export
MIData <- function(peak_areas, exp_names)
{
  # verify dimensions
  stopifnot(ncol(peak_areas) == length(exp_names) + 1)
  # create object
  mi_data <- list()
  class(mi_data) <-"MIData"

  # temp unique peak ids
  peak_ids <- unique(peak_areas[[1]])

  # find zero carbon peaks
  # TODO: this introduces indexing errors.
  # It doesn't seem like a good idea to mutate the peak_areas matrix,
  # should be moved elsewhere
#  nr_carbon <- as.numeric(table(factor(peak_areas[[1]], levels = unique(peak_ids)))) - 1
#  zero_carbon_peaks <- as.numeric(names(table(
#      factor(
#      peak_areas[[1]], levels = unique(peak_ids)))[which(nr_carbon == 0)]
#  ))
  # update peak_areas by removing the zero carbon peaks
#  peak_areas <- peak_areas[-which(peak_areas[[1]] %in% zero_carbon_peaks), ]

  # unique peak ids
  mi_data$peak_ids <- unique(peak_areas[[1]])

  # start index of each peak
  mi_data$peak_index <- match(mi_data$peak_ids, peak_areas[[1]])
  # find no. atoms for each peak (no. MIs = no.atoms + 1)
  mi_data$peak_n_atoms <-
    as.numeric(table(factor(peak_areas[[1]], levels = mi_data$peak_ids))) - 1
  # # precompute list of peak index vectors for each atom number
  mi_data$n_atoms_index <- create_atom_index(mi_data$peak_n_atoms)

  # names of the tracing experiments
  mi_data$experiments <- unique(exp_names)
  # index of first replicate for each experiment
  mi_data$exp_index <- match(mi_data$experiments, exp_names)
  # number of replicates per experiment
  mi_data$exp_n_rep <-
    as.numeric(table(factor(exp_names, levels = mi_data$experiments)))

  # store the peak area data as matrix
  mi_data$peak_areas <- as.matrix(peak_areas[2:ncol(peak_areas)])
  # ensure all peak areas are non-negative
  stopifnot(min(mi_data$peak_areas) >= 0)

  # MIDs
  mi_data$mids <- matrix(
    nrow = nrow(mi_data$peak_areas), ncol = ncol(mi_data$peak_areas))
  # averaged MIDs (collapsing across replicates)
  mi_data$avg_mids <- matrix(
    nrow = nrow(mi_data$peak_areas), ncol = length(mi_data$experiments))
  # compute MIDs
  for(p in 1:length(mi_data$peak_ids)) {
    rows <- get_mi_indices(mi_data, p)
    for(e in 1:length(mi_data$experiments)) {
      cols <- get_exp_indices(mi_data, e)
      # compute MID
      pa <- mi_data$peak_areas[rows, cols, drop = FALSE]
      # normalize nonzero mids
      mi_data$mids[rows, cols] <- normalize_mids(pa)

      # corresponding averaged MIDs
      #
      # DS: we currently have a condition for collapsing replicates (as follows):
      # 1. average only nonzero replicates
      #
      mi_data$avg_mids[rows, e] <- collapse_replicates(mi_data$mids[rows, cols, drop = F])
      # mi_data$avg_mids[rows, e] <- rowMeans(mi_data$mids[rows, cols, drop = FALSE])
    }
  }
  return(mi_data)
}


create_atom_index <- function(peak_n_atoms)
{
  index <- lapply(
    unique(peak_n_atoms),
    function(n) which(peak_n_atoms == n))
  names(index) <- as.character(unique(peak_n_atoms))
  return(index)
}

midata_subset <- function(midata, peak_index)
{
  new_midata <- midata
  # subset peaks ...
  
  # recompute the atoms index
  mi_data$n_atoms_index <- create_atom_index(mi_data$peak_n_atoms)
  
  # subset the peak area matrices ...
}

#
# normalize each column in a matrix of positive values
# so that each column sums to 1
# skip columns whose sum is zero
#
normalize_mids <- function(mids)
{
  normalized.mids <- matrix(0, nrow(mids), ncol(mids))
  for (i in 1:ncol(mids)){
    if (sum(mids[,i] != 0)){
      normalized.mids[,i] <- mids[,i]/sum(mids[,i])
    }
  }
  return(normalized.mids)
}

#
#+ average replicates that don't sum to zero, and only if 2/3 replicates don't sum to zero,
# if only one replicate does not sum to 0, return that replicate,
# else return zeros of size n+1 by 1, where n is the number of C atoms
#
collapse_replicates <- function(mid.matrix)
{
  col_sums <- colSums(mid.matrix)
  n_nonzero <- sum(col_sums > 0)
  return(rowSums(mid.matrix) / n_nonzero)
}

#
# get peak areas for a given peak p, experiment e
# (indices into the peak and experiment vectors)
#
get_peak_areas <- function(mi_data, p, e)
{
  rows <- get_mi_indices(mi_data, p)
  cols <- get_exp_indices(mi_data, e)
  return(as.matrix(
    mi_data$peak_areas[rows, cols], nrow = length(rows), ncol = length(col)))
}

#
# get MIDs, as above
#
get_mids <- function(mi_data, p, e)
{
  return(
    mi_data$mids[get_mi_indices(mi_data, p), get_exp_indices(mi_data, e), drop = FALSE])
}

#'Get an averaged MID vector
#'
#' for a given peak.
#'
#' @param mi_data an MIData object
#' @param p the peak index
#' @param e the experiment index
#' @export
get_avg_mid <- function(mi_data, p, e)
{
  return(mi_data$avg_mids[get_mi_indices(mi_data, p), e])
}

#' Get averaged MID vectors
#'
#' for a given peak, for all experiments
#'
#' @param mi_data an MIData object
#' @param p the peak index
#' @returns a matrix where each column is the MID from an experiment
#' @export
get_avg_mid_all <- function(mi_data, index)
{
  return(mi_data$avg_mids[get_mi_indices(mi_data, index), ])
}

get_mi_indices <- function(mi_data, p)
{
  return(mi_data$peak_index[[p]] + 0:mi_data$peak_n_atoms[[p]])
}

get_exp_indices <- function(mi_data, e)
{
  return(mi_data$exp_index[[e]] + 0:(mi_data$exp_n_rep[[e]]-1))
}

#' Get the index of a list of peak identifers in an MIData object
#' @param mi_data an MIData object
#' @param peak_ids a list of peak identifiers
#' @export
get_peak_index <- function(mi_data, peak_ids)
{
  return(match(peak_ids, mi_data$peak_ids))
}

#' Get the index of peaks with a specific number of atoms
#' @param mi_data an MIData object
#' @param n_atoms number of atoms in peaks of interest
#' @export
get_peak_index_n_atoms <- function(mi_data, n_atoms)
{
  return(mi_data$n_atoms_index[[as.character(n_atoms)]])
}

