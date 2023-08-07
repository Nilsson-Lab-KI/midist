#
# MID data functions
#

#' Construct an MIData (mass isotopomer data) object
#'
#' @param peak_areas a peak_area data.frame. The first column of peak_areas
#' must be peak identifiers, repeated for each MI of the same peak.
#' The second column holds formulas (but )
#' @param exp_names a optional list of experiment names for columns 2,3... in peak_areas.
#' These names will replace the column names for those columns.
#' @export
MIData <- function(peak_areas, exp_names = NULL)
{
  # verify first two column names
  stopifnot(colnames(peak_areas)[1:2] == c("Metabolite", "Formula"))

  # check for experiment names
  if(!is.null(exp_names)) {
    # verify length matches
    stopifnot(ncol(peak_areas) - 2 == length(exp_names))
  }
  else {
    # use data frame column names as experiment names
    exp_names <- colnames(peak_areas)[-(1:2)]
  }

  # create object
  mi_data <- list()
  class(mi_data) <- "MIData"

  # unique peak ids
  mi_data$peak_ids <- unique(peak_areas[["Metabolite"]])
  # start index of each peak in the mass isotopomer data
  mi_data$peak_index <- match(mi_data$peak_ids, peak_areas[[1]])
  # unique formulas
  mi_data$peak_formulas <- peak_areas[["Formula"]][mi_data$peak_index]
  # find no. atoms for each peak (no. MIs = no.atoms + 1)
  mi_data$peak_n_atoms <-
    as.numeric(table(factor(peak_areas[["Metabolite"]], levels = mi_data$peak_ids))) - 1
  # verify that peak_n_atoms matches peak_index
  stopifnot(all(mi_data$peak_index == find_mi_index(mi_data$peak_n_atoms)))
  # # precompute list of peak index vectors for each atom number
  mi_data$n_atoms_index <- create_atom_index(mi_data$peak_n_atoms)

  # names of the tracing experiments
  mi_data$experiments <- unique(exp_names)
  # index of first replicate for each experiment
  mi_data$exp_index <- match(mi_data$experiments, exp_names)
  # number of replicates per experiment
  mi_data$exp_n_rep <-
    as.numeric(table(factor(exp_names, levels = mi_data$experiments)))

  # remove first two columns (peak_ids / metabolite names, formulas)
  peak_areas <- as.matrix(peak_areas[3:ncol(peak_areas)])

  # MIDs
  mi_data$mids <- matrix(
    nrow = nrow(peak_areas), ncol = ncol(peak_areas)
  )
  # compute MIDs
  for (p in 1:length(mi_data$peak_ids)) {
    rows <- get_mi_indices(mi_data, p)
    for (e in 1:length(mi_data$experiments)) {
      cols <- get_exp_indices(mi_data, e)
      # compute MID
      pa <- peak_areas[rows, cols, drop = FALSE]
      # normalize nonzero mids
      mi_data$mids[rows, cols] <- normalize_mids(pa)
    }
  }
  # averaged MIDs
  mi_data$avg_mids <- calc_avg_mids(mi_data)

  return(mi_data)
}


#' Read peak areas from file and create an MIData object
#'
#' @param peak_areas_fname Name of a a tab-separated file containing peak areas
#' @returns An MIData object
#' @export
get_midata <- function(peak_areas_fname)
{
  peak_areas <- as.data.frame(
    utils::read.delim(peak_areas_fname, header = T, sep = "\t", check.names = F))

  if ("MassIsotopomer" %in% colnames(peak_areas)) {
    peak_areas <- peak_areas[, -which(colnames(peak_areas) == "MassIsotopomer")]
  }

  # make an MIData object from peak_areas, fetching experiments from column names of the matrix
  midata <- MIData(peak_areas)

  return(midata)
}


#' Compute an index vector into the MI data table for given of atom sizes
#' @param peak_n_atoms A vector of atom sizes
#' @returns A vector of indices to the first peak for each given atom size
find_mi_index <- function(peak_n_atoms)
{
  n <- length(peak_n_atoms)
  return(c(1, (cumsum(peak_n_atoms + 1) + 1)[-n]))
}


#' Create an index list mapping each number of atoms n to the indices of the peaks
#' having n atoms
#'
#' NOTE: this is a more generic function, applicable to any list,
#' could be moved to util.R
#'
#' @param peak_n_atoms number of atoms per peak
create_atom_index <- function(peak_n_atoms)
{
  index <- lapply(
    unique(peak_n_atoms),
    function(n) which(peak_n_atoms == n)
  )
  names(index) <- as.character(unique(peak_n_atoms))
  return(index)
}


# compute average MIDs (collapsing across replicates)
calc_avg_mids <- function(mi_data)
{
  avg_mids <- matrix(
    nrow = nrow(mi_data$mids), ncol = length(mi_data$experiments)
  )
  # compute MIDs
  for (p in 1:length(mi_data$peak_ids)) {
    rows <- get_mi_indices(mi_data, p)
    for (e in 1:length(mi_data$experiments)) {
      cols <- get_exp_indices(mi_data, e)
      # collapse across replicates
      avg_mids[rows, e] <- collapse_replicates(mi_data$mids[rows, cols, drop = F])
    }
  }
  return(avg_mids)
}


#' Subset an MIData object to the peaks given by peak_index, and return a new MIData object
#'
#' @param mi_data MIData object
#' @param peak_subset_index indices of peaks to be included in the subset
#' @export
midata_subset <- function(mi_data, peak_subset_index)
{
  # create subset midata object
  midata_subset <- list()
  class(midata_subset) <- "MIData"

  new_n_peaks <- length(peak_subset_index)
  new_n_atoms <- mi_data$peak_n_atoms[peak_subset_index]

  # NOTE: the order or assigning fields below matters,
  # as fields get assigned a position 1,2,... depending on order,
  # and e.g. testthat::expect_equal will fail if the positions differ

  # unique peak ids
  midata_subset$peak_ids <- mi_data$peak_ids[peak_subset_index]

  # index of each peak into the new MI data table
  midata_subset$peak_index <- find_mi_index(new_n_atoms)
  # unique peak formulas
  midata_subset$peak_formulas <- mi_data$peak_formulas[peak_subset_index]
  # number of atoms per peak
  midata_subset$peak_n_atoms <- new_n_atoms
  # indices of peaks per C group
  midata_subset$n_atoms_index <- create_atom_index(midata_subset$peak_n_atoms)

  # experiments are unchanged
  midata_subset$experiments <- mi_data$experiments
  midata_subset$exp_index <- mi_data$exp_index
  midata_subset$exp_n_rep <- mi_data$exp_n_rep

  # subset mids and avg_mids
  mi_index <- unlist(lapply(peak_subset_index,
                            function(x) get_mi_indices(mi_data, x)))
  midata_subset$mids <- mi_data$mids[mi_index, , drop = FALSE]
  midata_subset$avg_mids <- mi_data$avg_mids[mi_index, , drop = FALSE]

  return(midata_subset)
}


#' Apply a function to each MID in an MIData object
#'
#' This can be used for example to apply 13C correction to a data set
#' @param midata An MIData object
#' @param f The function to apply. Must take an MID vector as first argument
#' and return a valid MID vector of the same size.
#' @export
#'
midata_transform <- function(midata, f)
{
  # copy MIData object
  new_midata <- midata
  # apply function to each peak and sample
  n_col <- ncol(midata$mids)
  for (p in 1:length(midata$peak_ids)) {
    rows <- get_mi_indices(midata, p)
    for (i in 1:n_col) {
      new_midata$mids[rows, i] <- f(midata$mids[rows, i])
    }
  }
  # updata the averaged MIDs
  new_midata$avg_mids <- calc_avg_mids(new_midata)
  return(new_midata)
}


#' Reorder the peak IDs of an mi_data object
#'
#' @param midata An MIData object
#' @returns A vector of reordered peak IDs
#' @export
misplace_peak_ids <- function(midata)
{
  for (i in 1:length(midata$n_atoms_index)) {
    if (length(midata$n_atoms_index[[i]]) != 1) {
      # Shuffle the vector while ensuring none of the elements are in their original spot
      # Note that this is not a complete randomization, but rather misplacing every element
      new_ind <- midata$n_atoms_index[[i]]
      while (any(new_ind == midata$n_atoms_index[[i]])) {
        new_ind <- sample(midata$n_atoms_index[[i]])
      }
      # update the order of peak IDs without changing MIDs
      midata$peak_ids[midata$n_atoms_index[[i]]] <- midata$peak_ids[new_ind]
      rm(new_ind)
    } else midata$peak_ids[midata$n_atoms_index[[i]]] <- midata$peak_ids[midata$n_atoms_index[[i]]]

  }
  return(midata$peak_ids)
}


#' Censor false mass isotopomers
#'
#' Find mass isotopomers whose fraction is above the given threshold
#' after correction for natural 13C in at least the given number of
#' experiments (inclusive), set them to zero, and renormalize
#' @param mi_data An MIData object
#' @param threshold MI fraction threshold (after 13C correction)
#' @param min_experiments Minimum number of experiments where an MI is observed
#' to be considered false
#' @returns A new MIData object with false MIs censored
#' @export
censor_false_mi <- function(mi_data, threshold = 0.03, min_experiments = 1)
{
  # copy MIData object
  new_midata <- mi_data

  for (p in 1:length(mi_data$peak_ids)) {
    # get averaged mids for this peak, for all experiments
    avg_mids <- get_avg_mid(mi_data, p)

    # correct these mids for naturally occurring 13-C
    corrected_mids <- apply(avg_mids, 2, c13correct)

    # find index of false MIs
    false_index <- find_false_mi(corrected_mids, threshold, min_experiments)

    # rows of this peak
    rows <- get_mi_indices(mi_data, p)
    # replace false isotopes by zero
    new_mids <- new_midata$mids[rows, ]
    new_mids[false_index,] <- 0
    # renormalize
    new_midata$mids[rows, ] <- t(t(new_mids) / colSums(new_mids))
  }
  # update the averaged MIDs
  new_midata$avg_mids <- calc_avg_mids(new_midata)
  return(new_midata)
}

find_false_mi <- function(corrected_mids, threshold, min_experiments)
{
  # count number of experiments with an MI above threshold
  n_above_threshold <- as.vector(
    apply(corrected_mids[-1, , drop = FALSE], 1, function(x) sum(x > threshold)))

  # find indices of "false" isotopes and add 1 to account for M+0
  return(which(n_above_threshold >= min_experiments) + 1)
}


#
# normalize each column in a matrix of positive values
# so that each column sums to 1
# skip columns whose sum is zero
# TODO: move this to mid.R ?
#
normalize_mids <- function(mids)
{
  normalized.mids <- matrix(0, nrow(mids), ncol(mids))
  for (i in 1:ncol(mids)) {
    if (sum(mids[, i]) != 0 & unique(is.na(mids[,i])) == F) {
      normalized.mids[, i] <- mids[, i] / sum(mids[, i])
    }
  }
  return(normalized.mids)
}


#
# average replicates, excluding any MIDs that sum to zero
# else return zeros of size n+1 by 1, where n is the number of C atoms
#
collapse_replicates <- function(mid_matrix)
{
  col_sums <- colSums(mid_matrix)
  n_nonzero <- sum(col_sums > 0)
  if (n_nonzero > 0) {
    return(rowSums(mid_matrix) / n_nonzero)
  } else {
    return(rep(0, nrow(mid_matrix)))
  }
}


#' Get MIDs for a given peak and experiment from an MIData object
#'
#' @param mi_data an MIData object
#' @param p the peak index
#' @param e the experiment index
#' @returns A matrix with MIDs in columns
#' @export
get_mids <- function(mi_data, p, e)
{
  return(
    mi_data$mids[get_mi_indices(mi_data, p), get_exp_indices(mi_data, e), drop = FALSE]
  )
}


#' Get averaged MIDs from an MIData object
#'
#' Retrieve a vector of averages MIDs for the given peak p and experiment e;
#' or, if e is omitted, the matrix of MIDs for peak p across all experiments.
#'
#' @param mi_data an MIData object
#' @param p the peak index
#' @param e the experiment index (optional)
#' @return An MID vector, or, if e is omitted, a matrix whose columns are MIDs
#' @export
get_avg_mid <- function(mi_data, p, e)
{
  return(mi_data$avg_mids[get_mi_indices(mi_data, p), e])
}


#' Get MI indices of a given peak
#'
#' @param mi_data an MIData object
#' @param p the peak index
#' @returns a vector of MI indices
get_mi_indices <- function(mi_data, p)
  {
  return(mi_data$peak_index[[p]] + 0:mi_data$peak_n_atoms[[p]])
}


#' Get indices into the columns of the $mid matrix for experiment e
#' @param mi_data An MIData object
#' @param e An experiment index
#' @export
get_exp_indices <- function(mi_data, e) {
  return(mi_data$exp_index[[e]] + 0:(mi_data$exp_n_rep[[e]] - 1))
}


#' Get the index of a list of peak identifers in an MIData object
#' @param mi_data an MIData object
#' @param peak_ids a list of peak identifiers
#' @export
get_peak_index <- function(mi_data, peak_ids) {
  return(match(peak_ids, mi_data$peak_ids))
}


#' Get the index of peaks with a specific number of atoms
#' @param mi_data an MIData object
#' @param n_atoms number of atoms in peaks of interest
get_peak_index_n_atoms <- function(mi_data, n_atoms) {
  return(mi_data$n_atoms_index[[as.character(n_atoms)]])
}


#' Get the number of atoms for a given peak
#' @param mi_data an MIData object
#' @param p the peak index
#' @export
get_peak_n_atoms <- function(mi_data, p) {
  return(mi_data$peak_n_atoms[[p]])
}

