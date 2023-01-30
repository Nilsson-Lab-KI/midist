#
# MID data functions
#

#' Construct an MIData (mass isotopomer data) object
#'
#' @param peak_areas a peak_area data.frame. The first column of peak_areas must be peak identifiers,
#' repeated for each MI of the same peak, and MIs must be increasing 0,1,...,n for each peak
#' @param exp_names a optional list of experiment names matching columns 2,3... in peak_areas
#' @export
MIData <- function(peak_areas, exp_names) {
  # # verify dimensions
  # stopifnot(ncol(peak_areas) == length(exp_names) + 1)

  # verify first three column names
  stopifnot(colnames(peak_areas)[1:2] %in% c("Metabolite", "Formula"))

  # create object
  mi_data <- list()
  class(mi_data) <- "MIData"

  # unique peak ids
  mi_data$peak_ids <- unique(peak_areas[["Metabolite"]])
  # start index of each peak
  mi_data$peak_index <- match(mi_data$peak_ids, peak_areas[[1]])
  # unique formulas
  mi_data$peak_formulas <- peak_areas[["Formula"]][mi_data$peak_index]
  # find no. atoms for each peak (no. MIs = no.atoms + 1)
  mi_data$peak_n_atoms <-
    as.numeric(table(factor(peak_areas[["Metabolite"]], levels = mi_data$peak_ids))) - 1
  # # precompute list of peak index vectors for each atom number
  mi_data$n_atoms_index <- create_atom_index(mi_data$peak_n_atoms)

  # names of the tracing experiments
  mi_data$experiments <- unique(exp_names)
  # index of first replicate for each experiment
  mi_data$exp_index <- match(mi_data$experiments, exp_names)
  # number of replicates per experiment
  mi_data$exp_n_rep <-
    as.numeric(table(factor(exp_names, levels = mi_data$experiments)))

  # remove first three columns (peak_ids / metabolite names, formulas, and mass isotopomers)
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


#' Create an index list mapping each number of atoms n to the indices of the peaks having n atoms
#'
#' @param peak_n_atoms number of C atoms per peak
#' @export
create_atom_index <- function(peak_n_atoms) {
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
      avg_mids[rows, e] <- collapse_replicates(mi_data$mids[rows, cols, drop = F])
    }
  }
  return(avg_mids)
}


#
#
# TODO - WORK MORE ON THIS

#' Subset an MIData object to the peaks given by peak_index, and return a new MIData object
#'
#' @param mi_data MIData object
#' @param peak_index indices of peaks to be included in the subset
#' @export
midata_subset <- function(mi_data, peak_index) {
  # create subset midata object
  midata_subset <- list()
  class(midata_subset) <- "MIData"

  # unique peak ids
  midata_subset$peak_ids <- mi_data$peak_ids[peak_index]

  # unique peak formulas
  midata_subset$peak_formulas <- mi_data$peak_formulas[peak_index]

  # number of atoms per peak
  midata_subset$peak_n_atoms <- mi_data$peak_n_atoms[peak_index]

  # indices of peaks per C group
  midata_subset$n_atoms_index <- create_atom_index(midata_subset$peak_n_atoms)

  # peak indices
  midata_subset$peak_index <- match(midata_subset$peak_ids, rep(midata_subset$peak_ids, (midata_subset$peak_n_atoms + 1)))

  # experiments
  midata_subset$experiments <- mi_data$experiments

  # experiment index
  midata_subset$exp_index <- match(midata_subset$experiments, as.vector(unique(as.factor(midata_subset$experiments))))

  midata_subset$exp_n_rep <-
    as.numeric(table(factor(midata_subset$experiments, levels = midata_subset$experiments)))

  # subset mids and avg_mids
  midata_subset$mids <- as.matrix(do.call(rbind.data.frame, lapply(lapply(peak_index, function(x, mi_data) get_avg_mid(mi_data, x), mi_data), function(x) as.matrix(x))))
  colnames(midata_subset$mids) <- rownames(midata_subset$mids) <- NULL
  midata_subset$avg_mids <- as.matrix(do.call(rbind.data.frame, lapply(lapply(peak_index, function(x, mi_data) get_avg_mid(mi_data, x), mi_data), function(x) as.matrix(x))))
  colnames(midata_subset$avg_mids) <- rownames(midata_subset$avg_mids) <- NULL

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

#
# normalize each column in a matrix of positive values
# so that each column sums to 1
# skip columns whose sum is zero
#
normalize_mids <- function(mids) {
  normalized.mids <- matrix(0, nrow(mids), ncol(mids))
  for (i in 1:ncol(mids)) {
    if (sum(mids[, i] != 0)) {
      normalized.mids[, i] <- mids[, i] / sum(mids[, i])
    }
  }
  return(normalized.mids)
}

#
# average replicates, excluding any MIDs that sum to zero
# else return zeros of size n+1 by 1, where n is the number of C atoms
#
collapse_replicates <- function(mid_matrix) {
  col_sums <- colSums(mid_matrix)
  n_nonzero <- sum(col_sums > 0)
  if (n_nonzero > 0) {
    return(rowSums(mid_matrix) / n_nonzero)
  } else {
    return(rep(0, nrow(mid_matrix)))
  }
}


#
# get MIDs, as above
#
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
#' for a given peak, for all experiments
#'
#' @param mi_data an MIData object
#' @param index the peak index
#' @returns a matrix where each column is the MID from an experiment
#' @export
get_avg_mid_all <- function(mi_data, index) {
  return(mi_data$avg_mids[get_mi_indices(mi_data, index), ])
}


#' @export
get_avg_mids_by_size <- function(mi_data, n_atoms, e) {
  index <- get_peak_index_n_atoms(mi_data, n_atoms)
  return(sapply(index, function(i) get_avg_mid(mi_data, i, e)))
}


#' Get MI indices of a given peak
#'
#' @param mi_data an MIData object
#' @param p the peak index
#' @returns a vector of MI indices
#' @export
get_mi_indices <- function(mi_data, p) {
  return(mi_data$peak_index[[p]] + 0:mi_data$peak_n_atoms[[p]])
}

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
#' @export
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

#' @export
get_formula <- function(mi_data, p) {
  return(
    mi_data$peak_formulas[p]
  )
}

#' @export
get_mass <- function(midata, p) {
  return(
    midata$peak_masses[p]
  )
}
