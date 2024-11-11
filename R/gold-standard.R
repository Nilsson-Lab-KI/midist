#
# gold standard functions
#


#' Compute accuracy measures element-wise for two binary matrices
#' @param true_matrix Binary matrix where 1 indicates true
#' @param predicted_matrix Binary matrix where 1 indicates a positive prediction
#' @returns A list with elements $tpr, $fpr, $precision, $fdr, $f1_score
calc_accuracy <- function(true_matrix, predicted_matrix)
{
  accuracy <- list()
  # basic counts
  tp <- length(which(true_matrix == 1 & predicted_matrix == 1))
  fp <- length(which(true_matrix == 0 & predicted_matrix == 1))
  tn <- length(which(true_matrix == 0 & predicted_matrix == 0))
  fn <- length(which(true_matrix == 1 & predicted_matrix == 0))
  stopifnot(tp+fp+tn+fn == nrow(true_matrix)*ncol(true_matrix))
  #
  # accuracy measures
  #
  accuracy$tpr <- tp/(tp+fn)
  accuracy$fpr <- fp/(fp+tn)
  accuracy$precision <- tp/(tp+fp)
  accuracy$fdr <- fp/(fp+tp)
  accuracy$f1_score <- 2*tp/(2*tp + fp + fn)

  return(accuracy)
}

#' Compute accuracy measures at a given percentile cutoffs for pairwise_matrix

#' @param pairwise_matrix A symmetric matrix of distances
#' @param gold_standard A symmetric, binary matrix where 1 indicates a true pair
#' @param measure ?
#' @param noise ?
#' @param experiment ?
#' @param subset_size ?
#' @param percentiles A vector of percentiles, in the range 0--100
#'
#' #' @returns A data frame with accuracy values at each percentile
#' @export
get_global_percentile_accuracy <- function(
    pairwise_matrix, gold_standard, measure, noise, experiment, subset_size, percentiles)
{
  # create an empty data frame to store the accuracy results
  accuracy_df <- data.frame(measure = NA, noise = NA, percentile = NA,
                            tpr = NA, fpr = NA, precision = NA, fdr = NA, f1_score = NA,
                            experiment = NA, subset_size = NA)

  # make sure the diagonal is NA
  diag(pairwise_matrix) <- NA

  # now we go from full to empty
  tpr <- c(); fpr <- c(); precision <- c(); fdr <- c(); f1_score <- c()
  for (p in 1:length(percentiles)){
    # filter for percentile
    filtered_pm <- filter_pairwise_matrix_global(pairwise_matrix, percentile = percentiles[p]*0.01)
    # binarize
    filtered_pm[which(is.na(filtered_pm))] <- 0
    filtered_pm[which(filtered_pm != 0)] <- 1

    # now we have two binary matrices to compare to each other: from gs-fraction and MID-distance
    accuracy <- calc_accuracy(gold_standard, filtered_pm)
    tpr[p] <- accuracy$tpr
    fpr[p] <- accuracy$fpr
    precision[p] <- accuracy$precision
    fdr[p] <- accuracy$fdr
    f1_score[p] <- accuracy$f1_score
  }

  accuracy_df <- rbind(accuracy_df,
                       data.frame(
                         measure = measure, noise = noise, percentile = percentiles,
                         tpr = tpr, fpr = fpr, precision = precision, fdr = fdr, f1_score = f1_score,
                         experiment = experiment, subset_size = subset_size))

  # remove the NA line
  accuracy_df <- accuracy_df[-1,]

  return(accuracy_df)

}


#' Compute precision and recall for all pairs in a pairwise_matrix
#'
#' @param distance_matrix A symmetric matrix of distances
#' @param gold_standard A symmetric, binary matrix where 1 indicates a true pair
#' @returns A data.frame with columns 'precision' and 'recall'
#' @export
#'
continuous_accuracy <- function(distance_matrix, gold_standard)
{
  # vectors of lower-diagonal elements
  gs_elements <- gold_standard %>%
    get_lower_triangular()
  dm_elements <- distance_matrix %>%
    get_lower_triangular() %>%
    replace_na_with_max()

  # order by increasing distance
  dm_order_index <- order(dm_elements)
  gs_ordered <- gs_elements[dm_order_index]

  # cumulative sum of true and false pairs from gold standard
  tp <- cumsum(gs_ordered)
  fp <- cumsum(abs(1 - gs_ordered))

  # precision and recall for a single decision threshold
  return(
    data.frame(
      precision = tp / (tp + fp),
      recall = tp / max(tp)
    )
  )
}


#' Compute precision and recall for all pairs in a distance matrix
#'
#' @param distance_matrix A symmetric matrix of distances
#' @param gold_standard A symmetric, binary matrix where 1 indicates a true pair
#' @param measure ?
#' @param subset_size ?
#' @param subset_sample_no ?
#' @param noise ?
#' @param experiment ?
#' @export
get_continuous_accuracy <- function(distance_matrix, gold_standard,
                                    measure, subset_size, subset_sample_no,
                                    noise, experiment)
{
  # make sure the diagonal is NA
  diag(distance_matrix) <- NA
  diag(gold_standard) <- NA

  # put both matrices in a vector
  gs_vector <- as.vector(gold_standard)
  pm_vector <- as.vector(distance_matrix)

  # diagonals in gold standard are NAs - remove them from both vectors
  gs_na_ind <- which(is.na(gs_vector))
  pm <- pm_vector[-gs_na_ind]
  gs <- gs_vector[-gs_na_ind]

  # now assign a maximal distance to all remaining NAs of pm
  pm[which(is.na(pm))] <- max(pm, na.rm = T)

  # sort matrices by increasing distance
  order_ind <- order(pm)
  pm <- pm[order_ind]
  gs <- gs[order_ind]

  # cumulative sum of true and false pairs from gold standard
  tp <- cumsum(gs)
  fp <- cumsum(abs(1-gs))

  # Compute precision and recall for a single decision threshold
  threshold <- 1:length(gs)
  precision <- tp[threshold] / (tp[threshold] + fp[threshold])
  recall <- tp[threshold] / max(tp)


  return(data.frame(measure = measure,
                    subset_size = subset_size, subset_sample_no = subset_sample_no,
                    noise = noise,
                    tpr = recall, precision = precision,
                    experiment = experiment)
  )

}



