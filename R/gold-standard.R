# 
# gold standard functions
#
#' @export
read_gold_standard <- function(gs_fname, symmetric = T, binary = T, sparsity_cutoff = 0.5){
  
  gs <- as.matrix(read.delim(gs_fname, header = T, sep = "\t", check.names = F))
  
  if (symmetric == T){
    # symmetrize gs
    for (r2 in 1:nrow(gs)){
      for (d2 in r2:ncol(gs)){
        gs[r2,d2] <- gs[d2,r2] <- max(gs[r2,d2], gs[d2,r2])
      }
    }
  }
  
  # threshold gs based on the sparsity cutoff
  gs[which(gs < sparsity_cutoff)] <- 0
  
  if (binary == T)
    gs[which(gs != 0)] <- 1
  
  return(gs)
}

#' @export
# function to convolute: x + z -> y
convolute_gs <- function(z, x, y, gs){
  return((1-gs[z,x])*gs[x,y] + (1-gs[x,z])*gs[z,y])
}

#' @export
get_convoluted_gs <- function(mmm, gs_raw, input, symmetrize_by = max){
  # compare this combined and filtered pairwise matrix to the fractional gold standard
  unique_z <- unique(mmm[which(is.na(mmm) == F)])
  
  # convolute the gold standard for UNWEIGHTED
  gs <- gs_raw
  for (r2 in 1:nrow(gs_raw)){
    for (d2 in 1:nrow(gs_raw)){
      if (mmm[r2,d2] %in% unique_z){
        if (length(get_avg_mid(input$midata, r2, 1)) < length(get_avg_mid(input$midata, d2, 1))) 
          gs[r2,d2] <- convolute_gs(mmm[r2,d2], r2, d2, gs_raw) else
            gs[r2,d2] <- convolute_gs(mmm[r2,d2], d2, r2, gs_raw)
      }
    }
  }
  # now symmetric
  for (r2 in 1:nrow(gs)){
    for (d2 in r2:ncol(gs)){
      gs[r2,d2] <- gs[d2,r2] <- symmetrize_by(gs[r2,d2], gs[d2,r2])
    }
  }
  
  return(gs)
}

# some functions that I don't know where to put
#' @export
calc_accuracy <- function(true_matrix, predicted_matrix){
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
#' @export
get_percentile_accuracy <- function(pairwise_matrix, gold_standard, measure, noise, experiment, subset_size, percentiles){
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
    filtered_pm <- filter_pairwise_matrix(pairwise_matrix, percentile = percentiles[p]*0.01)
    # binarize
    filtered_pm[which(is.na(filtered_pm))] <- 0
    filtered_pm[which(filtered_pm != 0)] <- 1
    
    # now we have to binary matrices to compare to each other: from gs-fraction and MID-distance
    accuracy <- calc_accuracy(gold_standard, filtered_pm)
    tpr[p] <- accuracy$tpr
    fpr[p] <- accuracy$fpr
    precision[p] <- accuracy$precision
    fdr[p] <- accuracy$fdr
    f1_score[p] <- accuracy$f1_score
  }
  
  accuracy_df <- rbind(accuracy_df, data.frame(measure = measure, noise = noise, percentile = percentiles,
                                               tpr = tpr, fpr = fpr, precision = precision, fdr = fdr, f1_score = f1_score,
                                               experiment = experiment, subset_size = subset_size))
  
  # remove the NA line
  accuracy_df <- accuracy_df[-1,]
  
  return(accuracy_df)
  
}
#' @export
get_global_percentile_accuracy <- function(pairwise_matrix, gold_standard, measure, noise, experiment, subset_size, percentiles){
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
    
    # now we have to binary matrices to compare to each other: from gs-fraction and MID-distance
    accuracy <- calc_accuracy(gold_standard, filtered_pm)
    tpr[p] <- accuracy$tpr
    fpr[p] <- accuracy$fpr
    precision[p] <- accuracy$precision
    fdr[p] <- accuracy$fdr
    f1_score[p] <- accuracy$f1_score
  }
  
  accuracy_df <- rbind(accuracy_df, data.frame(measure = measure, noise = noise, percentile = percentiles,
                                               tpr = tpr, fpr = fpr, precision = precision, fdr = fdr, f1_score = f1_score,
                                               experiment = experiment, subset_size = subset_size))
  
  # remove the NA line
  accuracy_df <- accuracy_df[-1,]
  
  return(accuracy_df)
  
}



#' @export
get_continuous_accuracy <- function(pairwise_matrix, gold_standard, measure, noise, experiment, subset_size){

  # make sure the diagonal is NA
  diag(pairwise_matrix) <- NA
  diag(gold_standard) <- NA
  
  # put both matrices in a vector
  gs_vector <- as.vector(gold_standard)
  pm_vector <- as.vector(pairwise_matrix)
  
  # diagonals in gold standard are NAs - remove them from both vectors
  gs_na_ind <- which(is.na(gs_vector))
  pm <- pm_vector[-gs_na_ind]
  gs <- gs_vector[-gs_na_ind]
  
  # now assign a maximal distance to all remaining NAs of pm
  pm[which(is.na(pm))] <- max(pm, na.rm = T)
  
  # sort matrices
  order_ind <- order(pm)
  pm <- pm[order_ind]
  gs <- gs[order_ind]
  
  
  tp <- cumsum(gs)
  fp <- cumsum(abs(1-gs))
  # tn <- c(0:(length(gs)-1))*fp
  # fn <- c(0:(length(gs)-1))*tp
  
  # Compute precision and recall for a single decision threshold
  threshold <- 1:length(gs)
  precision <- tp[threshold] / (tp[threshold] + fp[threshold])
  recall <- tp[threshold] / max(tp)
  # precision <-  sapply(1:length(gs), function(x, tp, fp) tp[x] / (tp[x] + fp[x]), tp, fp) 
  # recall <-  sapply(1:length(gs), function(x, tp) tp[x] / max(tp), tp) 

  
  return(data.frame(measure = measure, noise = noise,
                    tpr = recall, precision = precision,
                    experiment = experiment, subset_size = subset_size))
  
}



#' @export
get_threshold_accuracy <- function(pairwise_matrix, gold_standard, measure, noise, experiment, thresholds){
  # create an empty data frame to store the accuracy results
  accuracy_df <- data.frame(measure = NA, noise = NA, threshold = NA, 
                            tpr = NA, fpr = NA, precision = NA, fdr = NA, f1_score = NA,
                            experiment = NA)
  
  # make sure the diagonal is NA
  diag(pairwise_matrix) <- NA
  pairwise_matrix[which(is.na(pairwise_matrix))] <- max(pairwise_matrix, na.rm = T)
  
  # now we go from full to empty
  tpr <- c(); fpr <- c(); precision <- c(); fdr <- c(); f1_score <- c()
  for (p in 1:length(thresholds)){
    # filter for threshold
    filtered_pm <- pairwise_matrix
    filtered_pm[which(filtered_pm > thresholds[p])] <- 0
    # binarize
    filtered_pm[which(filtered_pm != 0)] <- 1
    
    # now we have to binary matrices to compare to each other: from gs-fraction and MID-distance
    accuracy <- calc_accuracy(gold_standard, filtered_pm)
    tpr[p] <- accuracy$tpr
    fpr[p] <- accuracy$fpr
    precision[p] <- accuracy$precision
    fdr[p] <- accuracy$fdr
    f1_score[p] <- accuracy$f1_score
  }
  
  accuracy_df <- rbind(accuracy_df, data.frame(measure = measure, noise = noise, threshold = thresholds,
                                               tpr = tpr, fpr = fpr, precision = precision, fdr = fdr, f1_score = f1_score,
                                               experiment = experiment))
  
  # remove the NA line
  accuracy_df <- accuracy_df[-1,]
  
  return(accuracy_df)
  
}

  



