# LIBRARIES
library(lsei)


# CARBON NUMBER #
{
  # function for calculating the number of Carbon atoms for a given metabolite from MID data
  # requires a specific data format: a 2-column data frame whose column names are "Metabolite" and "MID"
  get_c_number <- function(metabolite_vector, mid_data){
    # make sure the MID data is useful
    stopifnot(length(which(is.na(match(c("Metabolite", "MID"), colnames(mid_data))))) == 0)
    
    mid_data <- mid_data[which(mid_data$Metabolite %in% metabolite_vector),]
    
    # factorising metabolite names so that the order does not get screwed up when counting the number of atoms
    mid_data$Metabolite <- factor(mid_data$Metabolite, levels = unique(mid_data$Metabolite))
    metabolite_vector <- factor(metabolite_vector, levels = unique(metabolite_vector))
    
    # counting the occurrence of each element subtracting 1 from which will yield the number of C atoms
    metabolites <- as.vector(as.data.frame(table(mid_data$Metabolite))$Var1)
    counts <- as.vector(as.data.frame(table(mid_data$Metabolite))$Freq)
    
    # finding metabolites in metabolite_vector, and fixing the order
    index <- match(metabolites, metabolite_vector)
    return(counts[index]-1)
    
  }
  
  # function for getting a matrix of carbon differences
  get_carbon_difference_matrix <- function(met_names, mid_data){
    
    # creating an empty matrix to be filled in with similarities/distances
    c_diff <- matrix(NA, length(met_names), length(met_names))
    colnames(c_diff) <- met_names
    rownames(c_diff) <- met_names
    
    for (r2 in 1:nrow(c_diff)){
      mid_r2 <- mid_data$MID[which(mid_data$Metabolite == met_names[r2])]
      nrcarbon_r2 <- length(mid_r2) - 1
      
      for (d2 in r2:ncol(c_diff)){
        print(c(r2,d2))
        mid_d2 <- mid_data$MID[which(mid_data$Metabolite == met_names[d2])]
        nrcarbon_d2 <- length(mid_d2) - 1
        
        c_diff[r2,d2] <- abs(nrcarbon_r2 - nrcarbon_d2)
        c_diff[d2,r2] <- abs(nrcarbon_r2 - nrcarbon_d2)
        
      }
      
    }
    
    diag(c_diff) <- 0
    return(c_diff)
    
  }
}


# ISOTOPIC ENRICHMENT #
{
  # calculate isotopic enrichment
  isotopic_enrichment <- function(mid) {
    n <- length(mid) - 1
    return(sum(mid * c(0:n)) / n)
  }
  
  # function for calculating isotopic enrichment of a metabolite for a particular experiment
  get_tracer_enrichment <- function(exp_name, met_subset){
    return(isotopic_enrichment(met_subset[, which(colnames(met_subset) == exp_name)]))
  }
  
  # function for getting enrichment profiles across tracers of each metabolite
  get_enrichment_profile <- function(met_name, sim_data, exp_names){
    met_subset <- sim_data[which(sim_data$Metabolite == met_name),]
    return(unlist(lapply(exp_names, get_tracer_enrichment, met_subset)))
  }
}


# CORRECTING NATURALLY OCCURRING 13-C #
{
  # this is a function to calculate the binomial distribution of a set of "i"s (B_i^n)
  # "is" is a vector of i values- that is [0, 1, 2, 3] for n = 3
  # p is the probability for the naturally occurring carbon-13, which is 0.0107 by default
  binomvals <- function(is, n, p){
    bins <- c()
    for (j in 1:length(is)){
      bin <- dbinom(is[j], n, p)
      bins <- c(bins, bin)
    }
    return(bins)
  }
  
  # correcting for the natural 13-C
  c13correct <- function(mid, p = 0.0107){
    
    nrCarbon <- length(mid) - 1
    end <- length(mid)
    correct <- matrix(0, end, end)
    
    for (d2 in 1:end){
      b1 <- c(0:(end-d2))
      b2 <- end-d2
      correct[d2:end, d2] <- binomvals(b1, b2, p)
    }
    
    return(pnnls(a = correct, b = mid)$x)
  }
  
}


# MID AND CONVOLUTION #
{
  bring_mid <- function(met, mid_data){
    return(mid_data$MID[which(mid_data$Metabolite == met)])
  }
  get_convoluted_mid <- function(middle_mid, s_mid){
    
    if (length(middle_mid) == length(s_mid) | length(middle_mid) > length(s_mid)){
      # get convolution matrix
      conv_matrix <- convolution_matrix_v2(middle_mid, s_mid)
      # convolute 
      c_mid <- conv_matrix %*% s_mid
      
      return(c_mid)
    }
    else {
      # get convolution matrix
      conv_matrix <- convolution_matrix_v2(s_mid, middle_mid)
      # convolute 
      c_mid <- conv_matrix %*% middle_mid
      
      return(c_mid)
    }
    
  }
  convolution_matrix_v2 <- function(longer_mid, shorter_mid){
    n_row <- length(longer_mid) + length(shorter_mid) - 1
    n_col <- length(shorter_mid)
    A <- matrix(0, n_row, n_col)
    
    for (i in 1:ncol(A)){
      A[i:(i+length(longer_mid)-1), i] <- longer_mid
    }
    
    return(A)
  }
}


# SIMILARITY #
{
  # cosine similarity
  cosine <- function(x, y){
    # cosi <- as.numeric(as.vector(unlist(x %*% y))) / (sqrt(sum(x^2)) * sqrt(sum(y^2)))
    cosi <- sum(x*y) / (sqrt(sum(x^2)) * sqrt(sum(y^2)))
    return(cosi)
  }
  
  # euclidean distance
  euclidean <- function(x, y) sqrt(sum((x - y)^2))
}


# DATA ANALYSIS #
{
  
}


# SIMILARITY/DISTANCE MATRIX ANALYSIS #
{
  
}


# PLOTTING #
{
  
}


