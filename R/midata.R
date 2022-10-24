#
# MID data functions
#

#
# Construct an MI (mass isotopomer) data object from a peak_area data.frame
# The first column of peak_areas must be peak identifiers, repeated for each
# MI of the same peak, and MIs must be increasing 0,1,...,n for each peak
# exp_name is a list of experiment names matching columns 2,3... in peak_areas
#
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
  nr_carbon <- as.numeric(table(factor(peak_areas[[1]], levels = unique(peak_ids)))) - 1
  zero_carbon_peaks <- as.numeric(names(table(factor(peak_areas[[1]], levels = unique(peak_ids)))[which(nr_carbon == 0)]))
  # update peak_areas by removing the zero carbon peaks
  peak_areas <- peak_areas[-which(peak_areas[[1]] %in% zero_carbon_peaks), ]

  # unique peak ids
  mi_data$peak_ids <- unique(peak_areas[[1]])

  # start index of each peak
  mi_data$peak_index <- match(mi_data$peak_ids, peak_areas[[1]])
  # find no. atoms for each peak (no. MIs = no.atoms + 1)
  mi_data$peak_n_atoms <-
    as.numeric(table(factor(peak_areas[[1]], levels = unique(mi_data$peak_ids)))) - 1

  # names of the tracing experiments
  mi_data$experiments <- unique(exp_names)
  # index of first replicate for each experiment
  mi_data$exp_index <- match(mi_data$experiments, exp_names)
  # number of replicates per experiment
  mi_data$exp_n_rep <-
    as.numeric(table(factor(exp_names, levels = unique(exp_names))))

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

#
# normalize mids, skip missing values to avoid division by zero
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
    mi_data$peak_areas[rows, cols], nrow = len(rows), ncol = len(col)))
}

#
# get MIDs, as above
#
get_mids <- function(mi_data, p, e)
{
  return(
    mi_data$mids[get_mi_indices(mi_data, p), get_exp_indices(mi_data, e), drop = FALSE])
}

#
# get an averaged MID vector for a given peak and experiment
#
get_avg_mid <- function(mi_data, p, e)
{
  return(mi_data$avg_mids[get_mi_indices(mi_data, p), e])
}

#
# get averaged MID vectors for a given peak, for all experiments
# returns an MI x experiments
get_avg_mid_all <- function(mi_data, p)
{
  return(mi_data$avg_mids[get_mi_indices(mi_data, p), ])
}

get_mi_indices <- function(mi_data, p)
{
  return(mi_data$peak_index[[p]] + 0:mi_data$peak_n_atoms[[p]])
}

get_exp_indices <- function(mi_data, e)
{
  return(mi_data$exp_index[[e]] + 0:(mi_data$exp_n_rep[[e]]-1))
}

get_peak_index <- function(mi_data, peak_ids)
{
  return(match(peak_ids, mi_data$peak_ids))
}

