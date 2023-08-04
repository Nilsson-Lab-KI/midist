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
      # find and replace zero peaks
      pa[, which(colSums(pa) == 0)] <- NA
      # normalize mids
      mi_data$mids[rows, cols] <- apply(pa, 2, function(vec) vec / sum(vec, na.rm = T))
    }
  }
  # averaged MIDs
  mi_data$avg_mids <- calc_avg_mids(mi_data)

  return(mi_data)
}


#' Read peak areas from file and create an MIData object
#'
#' @param peak_areas_fname Name of a a tab-separated file containing peak areas
#' @export
get_midata <- function(peak_areas_fname){

  peak_areas <- as.data.frame(read.delim(peak_areas_fname, header = T, sep = "\t", check.names = F))

  if ("MassIsotopomer" %in% colnames(peak_areas)) {
    peak_areas <- peak_areas[, -which(colnames(peak_areas) == "MassIsotopomer")]
  }

  # make an MIData object from peak_areas, fetching experiments from column names of the matrix
  midata <- MIData(peak_areas)

  return(midata)
}

#' Compute index vector into the MI data table
#' for a given list of atom sizes
find_mi_index <- function(peak_n_atoms)
{
  n <- length(peak_n_atoms)
  return(c(1, (cumsum(peak_n_atoms + 1) + 1)[-n]))
}



#' Create an index list mapping each number of atoms n to the indices of the peaks having n atoms
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
calc_avg_mids <- function(mi_data) {
  avg_mids <- matrix(
    nrow = nrow(mi_data$mids), ncol = length(mi_data$experiments)
  )
  # compute MIDs
  for (p in 1:length(mi_data$peak_ids)) {
    rows <- get_mi_indices(mi_data, p)
    for (e in 1:length(mi_data$experiments)) {
      cols <- get_exp_indices(mi_data, e)
      # collapse across replicates
      if (all(is.na(colSums(mi_data$mids[rows, cols, drop = F]))))
        avg_mids[rows, e] <- dbinom(c(0:(length(rows) - 1)), length(rows) - 1, 0.0107) else {
          new_areas <- rowSums(mi_data$mids[rows, cols, drop = F], na.rm = T) / length(which(!is.na(colSums(mi_data$mids[rows, cols, drop = F], na.rm = T))))
          # renormalization needed for some cases
          avg_mids[rows, e] <- new_areas / sum(new_areas)
        }
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

#'
#' Apply a function to each MID in an MIData object
#'
#' This can be used for example to apply 13C correction to a data set
#' @param mi_data An MIData object
#' @param f The function to apply. Must take an MID vector as first argument
#' and return a valid MID vector of the same size.
#' @export
#'
midata_transform <- function(midata, f) {
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

#'
#' Add replicates with some noise to an midata object where experiments are not repeated (and probably noise free also)
#' Returns a new midata object where the mids per peak per experiment are repeated with a defined level of noise
#' Updates all related fields in midata such as avg_mids, exp_index, etc.
#'
#' This is useful when the midata object contains simulated MIDs without replicates and noise,
#' and to add replicates we introduce a level of noise defined by the stdev parameter.
#'
#' @param midata An MIData object
#' @param stdev The desired standard deviation to generate dirichlet distribution
#' (see the random_mid() function for more details on the dirichlet distribution)
#' and return a valid MID vector of the same size.
#' @param nr_replicate Number of noisy replicates to generate
#'
add_noisy_replicates <- function(midata, stdev, nr_replicate) {
  # copy MIData object
  new_midata <- midata
  exp_names <- rep(midata$experiments, each = nr_replicate)

  # apply function to each peak
  new_midata$mids <- do.call(rbind.data.frame,
                             lapply(1:length(new_midata$peak_ids),
         function(p) do.call(cbind.data.frame, lapply(1:length(new_midata$experiments),
                                                      function(e) random_mid(get_avg_mid(new_midata, p, e),
                                                                             stdev, nr_replicate)))) )

  # need to replicate experiment names
  new_midata$exp_index <- match(new_midata$experiments, exp_names)
  # number of replicates per experiment
  new_midata$exp_n_rep <-
    as.numeric(table(factor(exp_names, levels = new_midata$experiments)))

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


midata_randomize <- function(midata){
  for (i in 1:length(midata$n_atoms_index)){
    print(i)
    if (length(midata$n_atoms_index[[i]]) != 1){
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
  return(midata)
}


#' Censor false mass isotopomers
#'
#' Find mass isotopomers whose fraction is above the given threshold
#' after correction for natural 13C in at least the given number of
#' experiments (inclusive), set them to zero, and renormalize
#' @param mi_data An MIData object
#' @param threshold MI fraction threshold (after 13C correction)
#' @param min_experiment Minimum number of experiments where an MI is observed
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
    
    # # NaN filtering here
    # nan <- all(is.na(colSums(new_midata$mids[rows, ])))
    # if (nan == F){
    #   nan_index <- which(is.na(colSums(new_midata$mids[rows, ])))
    #   new_midata$mids[rows, nan_index] <- dbinom(c(0:(length(rows) - 1)), length(rows) - 1, 0.0107)
    # }
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
# get MIDs, as above
#
#' @export
get_mids <- function(mi_data, p, e) {
  return(
    mi_data$mids[get_mi_indices(mi_data, p), get_exp_indices(mi_data, e), drop = FALSE]
  )
}


#' Get an averaged MID vector
#'
#' for a given peak.
#'
#' @param mi_data an MIData object
#' @param p the peak index
#' @param e the experiment index
#' @export
get_avg_mid <- function(mi_data, p, e) {
  return(mi_data$avg_mids[get_mi_indices(mi_data, p), e])
}

#' Get averaged MID vectors
#'
#' for a given peak index, for all experiments
#'
#' @param mi_data an MIData object
#' @param index the peak index
#' @returns a matrix where each column is the MID from an experiment
get_avg_mid_all <- function(mi_data, index) {
  return(mi_data$avg_mids[get_mi_indices(mi_data, index), ])
}


#' Get averaged MIDs for all metabolites of the given size
#' If no metabolites exists returns an empty list
#'
get_avg_mids_by_size <- function(mi_data, n_atoms, e) {
  index <- get_peak_index_n_atoms(mi_data, n_atoms)
  return(sapply(index, function(i) get_avg_mid(mi_data, i, e)))
}


#' Get MI indices of a given peak
#'
#' @param mi_data an MIData object
#' @param p the peak index
#' @returns a vector of MI indices
get_mi_indices <- function(mi_data, p) {
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


get_formula <- function(mi_data, p) {
  return(
    mi_data$peak_formulas[p]
  )
}


get_mass <- function(midata, p) {
  return(
    midata$peak_masses[p]
  )
}
