
### peak intensity processing functions:


calc_enrichment <- function(metabolite, tracer, mid_data, enrichment_tol = 0.0107)
{
  mid <- mid_data[which(mid_data$Metabolite == metabolite), which(colnames(mid_data) == tracer)]
  df <- data.frame(metabolite = metabolite, enrichment = isotopic_enrichment(mid))
  return(df)
}

get_enriched_mets <- function(metabolite, tracer, mid_data, enrichment_tol = 0.0107)
{
  mid <- mid_data[which(mid_data$Metabolite == metabolite), which(colnames(mid_data) == tracer)]
  enrichment <- isotopic_enrichment(mid)
  if (enrichment > enrichment_tol){
    return(metabolite)
  }
}

prepare <- function(mid_1, mid_2){
  if (length(mid_1) == length(mid_2)){
    return(list(mid_1, mid_2))
  }
  else if (length(mid_1) > length(mid_2)){
    return(list(mid_1, find_convolution(mid_1, mid_2)))
  }
  else
    return(list(find_convolution(mid_2, mid_1), mid_2))
}

correct_before_prepare <- function(mid_1, mid_2){
  mid_1_c <- c13correct(mid_1)
  mid_2_c <- c13correct(mid_2)
  if (length(mid_1_c) == length(mid_2_c)){
    return(list(mid_1_c[-1], mid_2_c[-1]))
  }
  else if (length(mid_1_c) > length(mid_2_c)){
    find_convolutiond_mid <- find_convolution(mid_1_c, mid_2_c)
    return(list(mid_1_c[-1], convoluted_mid[-1]))
  }
  else {
    convoluted_mid <- find_convolution(mid_2_c, mid_1_c)
    return(list(convoluted_mid[-1], mid_2_c[-1]))
  }

}

correct_after_prepare <- function(mid_1, mid_2){

  if (length(mid_1) == length(mid_2)){
    mid_1_c <- c13correct(mid_1)
    mid_2_c <- c13correct(mid_2)
    return(list(mid_1_c[-1], mid_2_c[-1]))
  }
  else if (length(mid_1) > length(mid_2)){
    convoluted_mid <- find_convolution(mid_1, mid_2)
    mid_1_c <- c13correct(mid_1)
    convoluted_mid_c <- c13correct(convoluted_mid)
    return(list(mid_1_c[-1], convoluted_mid_c[-1]))
  }
  else {
    convoluted_mid <- find_convolution(mid_2, mid_1)
    mid_2_c <- c13correct(mid_2)
    convoluted_mid_c <- c13correct(convoluted_mid)
    return(list(convoluted_mid[-1], mid_2_c[-1]))
  }

}

# this is a function to calculate the binomial distribution of a set of "i"s (B_i^n)
# "is" is a vector of i values- that is [0, 1, 2, 3] for n = 3
# p is the probability for the naturally occurring carbon-13, which is commonly considered 0.01
binomvals <- function(is, n, p){
  bins <- c()
  for (j in 1:length(is)){
    bin <- dbinom(is[j], n, p)
    bins <- c(bins, bin)
  }
  return(bins)
}


#
# get the length of the MID vector (number of atoms + 1)
# for a given peak number, based on the peak.areas data
#
get.mid.length <- function(peak_nr, peak.areas)
{
  subset <- peak.areas[which(peak.areas$Peak.nr == as.numeric(peak_nr)), 3]
  return(length(subset))
}

#
# get the MID corresponding to the peak area vector of a particular
# replicate of a tracer. Sets any MID with enrichment < 0.01 to zero
#
get.replicate.mid <- function(replicate, tracer.subset)
{
  area <- tracer.subset[, as.numeric(replicate)]

  if (sum(area) != 0) {
    frac <- area/sum(area)
    mid <- filter_enrichment(frac)
  }
  else {
    mid <- rep(0, length(area))
  }
  return(mid)
}

#
# get MIDs for three replicates of one tracer,
# averaging replicates
#
get.tracer.mids <- function(tracer.ind, peak.subset, tracers)
{
  tracer.subset <- peak.subset[, which(colnames(peak.subset) == tracers[as.numeric(tracer.ind)])]
  mids <- lapply(1:3, get.replicate.mid, tracer.subset)
  sums <- unlist(lapply(mids, sum))

  nr.nonzero <- length(which(sums != 0))

  if (nr.nonzero > 1){
    mid <- rowMeans(cbind.data.frame(mids)[, which(sums != 0)])
  }
  else {
    mid <- rep(0, nrow(tracer.subset))
  }
  return(mid)
}

#
# get MID for all tracers
#
get.peak.mids <- function(peak_nr, peak.areas, tracers)
{
  peak.subset <- peak.areas[which(peak.areas$Peak.nr == as.numeric(peak_nr)), ]
  mids <- lapply(1:length(tracers), get.tracer.mids, peak.subset, tracers)
  return(mids)
}


### dataframe creation functions:

#
# Create a vector of peak pair strings of the format "nnn_nnn"
# Any peaks in the "all.peaks" list preceding "peak" will be excluded,
# e.g. get_df(3, c(1,2,3,4,5)) => "3_3", "3_4", "3_5"
#
# TODO: rename this function
#
get.df <- function(peak, all.peaks)
{
  other.peaks <- all.peaks[which(all.peaks == peak):length(all.peaks)]
  df <- data.frame(met.1 = rep(as.numeric(peak), length(other.peaks)), met.2 = as.numeric(other.peaks))
  set <- unlist(apply(df, 1, paste0, collapse = "_"))
  return(set)
}

### cosine functions




# DS: this function will explode in the future as we switch between different approaches to handle uninformative mids
# improvements and documentation are necessary
# always be aware of the zero assignment here!!!

# returns the convoluted mid for a particular tracer. longer.mids and shorter.mids are lists of mids for tracers
convolute.mid <- function(tracer.ind, longer.mids, shorter.mids)
{
  longer.mid <- longer.mids[[tracer.ind]]

  # the other peak would be the smaller carbon one
  shorter.mid <- shorter.mids[[tracer.ind]]

  # need to have this function very fast
  solution.mid <- solution(longer.mid, shorter.mid)

  solution.enrichment <- filter_enrichment(solution.mid)

  if (sum(solution.enrichment) == 0){
    convoluted.mid <- rep(0, length(longer.mid))
  }

  else {
    convolution <- find_convolution(longer.mid, shorter.mid)
    convolution.enrichment <- filter_enrichment(convolution)

    if (sum(convolution.enrichment) == 0){
      convoluted.mid <- rep(0, length(longer.mid))
    }
    else {
      convoluted.mid <- convolution
    }
  }
  return(as.vector(convoluted.mid))
}

# this function is almost the same as solution(), but it can handle multiple tracers
# returns the solution mid for a particular tracer. longer.mids and shorter.mids are lists of mids for tracers
solution.mid <- function(tracer.ind, longer.mids, shorter.mids)
{

  longer.mid <- longer.mids[[tracer.ind]]

  # the other peak would be the smaller carbon one
  shorter.mid <- shorter.mids[[tracer.ind]]

  # need to have this function very fast
  solution.mid <- solution(longer.mid, shorter.mid)
  return(as.vector(solution.mid))
}

# returns true if neither MID for the given tracer is the zero vector
check.tracer.overlap <- function(tracer.index, mid.lists)
{
  mid.1 <- mid.lists[[1]][[tracer.index]]
  mid.2 <- mid.lists[[2]][[tracer.index]]
  return (sum(mid.1) != 0 & sum(mid.2) != 0)
}

# Given a pair of MID lists over tracers, this unlists (flattens)
# the lists, including only those tracers where both MIDs are labelled
# The output elements are ready to be used as input to calc.cosine()
prepare.mids <- function(mid.lists)
{
  # number of carbons of each MID in the pair
  carbons <- c(length(mid.lists[[1]][[1]])-1, length(mid.lists[[2]][[1]])-1)
  if (carbons[1] - carbons[2] == 0) {
    # equal number of carbons
    inds <- unlist(lapply(1:length(mid.lists[[1]]), check.tracer.overlap, mid.lists))
    mids <- list(unlist(mid.lists[[1]][inds]), unlist(mid.lists[[2]][inds]))
  }
  else {
    # unequal number of carbons
    longer.mid <- mid.lists[[which.max(carbons)]]
    shorter.mid <- mid.lists[[which.min(carbons)]]

    convoluted.mid <- lapply(1:length(mid.lists[[1]]), convolute.mid, longer.mid, shorter.mid)

    inds <- unlist(lapply(1:length(mid.lists[[1]]), check.tracer.overlap, list(longer.mid, convoluted.mid)))

    mids <- list(unlist(longer.mid[inds]), unlist(convoluted.mid[inds]))
  }
  return(mids)
}


# the main function that takes a peak pair, prepares mids (including convolutions),
# and output cosine scores for the input peak pair
get.cosine <- function(peak_pair, peak.mids)
{
  # get pair of MID data corresponding to the peak pair
  spl <- unlist(str_split(peak_pair, "_"))
  mid.lists <- list(peak.mids[[as.numeric(spl[1])]], peak.mids[[as.numeric(spl[2])]])
  # calculate convolutions if necessary
  # and flatten mids
  mids <- prepare.mids(mid.lists)
  cosine <- cosine_sim(mids[[1]], mids[[2]])
  return(cosine)
}

#' Apply a distance function dist to MID x,y after convolution
#'
#' if length(x) < length(y), we find a vector v such that
#' ||x * v - y|| is minimal, and compute dist(x*v, y)
#' and vice versa if length(x) > length(y)
#'
#' @param dist a distance function
#' @param x an MID
#' @param y an MID
#' @param conv_fn a function used to convolute MIDs
#' @returns the distance value
#' @export
get_mid_dist <- function(dist, x, y, conv_fn)
{
  # if there are any N/A values, the distance is N/A
  if (any(is.na(x)) || any(is.na(y)))
    return(NA)
  nx <- length(x)
  ny <- length(y)
  if (nx == ny)
    return(dist(x, y))
  else if (nx > ny) {
    conv <- conv_fn(x, y)
    return(dist(x, conv))
  }
  else {
    conv <- conv_fn(y, x)
    return(dist(y, conv))
  }

}



### post-cosine (evaluation) functions


# collects information for a given peak pair:
# compound formula, name (if available), cosine score, etc
get.peak.info <- function(index, network, peaks, names, formulas)
{
  print(index)
  indices <- as.numeric(unlist(str_split(network[index,1], "_")))
  peak.1 <- peaks[as.numeric(indices[1])]
  peak.2 <- peaks[as.numeric(indices[2])]
  formula.1 <- formulas[indices[1]]
  formula.2 <- formulas[indices[2]]
  name.1 <- names[indices[1]]
  name.2 <- names[indices[2]]
  cosine <- as.numeric(network[index,2])
  df <- data.frame(peak.id.1 = peak.1, formula.1 = formula.1, metabolite.1 = name.1,
                   peak.i.2 = peak.2, formula.2 = formula.2, metabolite.2 = name.2,
                   cosine)
  return(df)
}

get.carbon.number <- function(formula)
{
  # get the periodic table elements (symbols)
  Symbol <- c("C", "H", "He", "Li", "Be", "B", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu",
              "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce",
              "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr",
              "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv",
              "Ts", "Og")

  # make a data frame for the input formula
  brokenFormula <- as.data.frame(Symbol)
  brokenFormula$ElementNumber <- 0

  # splitting the formula into the smallest pieces
  formula_split <- unlist(strsplit(formula, ""))

  # merging the lower case letters (if there is any) to their upper case letter
  grep_lower <- grepl("[a-z]", formula_split)
  l <- which(grep_lower == TRUE)
  if (length(l) > 0){
    formula_split[l-1] <- paste(formula_split[l-1], formula_split[l], sep = "")
    formula_split[l] <- NA
    formula_split <- formula_split[-which(is.na(formula_split)==TRUE)]
  }

  # find elements' location in the broken formula
  element_ind <- which(formula_split %in% Symbol)

  # now we need to assign "1" in case the last element doesn't have a number
  if (element_ind[length(element_ind)] == length(formula_split)){
    formula_split <- c(formula_split, "1")
  }

  # now we need to assign "1" to any element in the middle that doesn't have a number
  if (length(element_ind) > 1){
    for (i in 1:(length(element_ind)-1)){
      row_index <- which(brokenFormula$Symbol == formula_split[element_ind[i]])
      if (element_ind[i]+1 == element_ind[i+1]){
        brokenFormula$ElementNumber[row_index] <- 1
      }
      if (element_ind[i]+1 != element_ind[i+1]){
        brokenFormula$ElementNumber[row_index] <- as.numeric(paste(formula_split[(element_ind[i]+1):(element_ind[i+1]-1)], collapse = ""))
      }
    }
  }

  i <- length(element_ind)
  row_index <- which(brokenFormula$Symbol == formula_split[element_ind[i]])
  brokenFormula$ElementNumber[row_index] <- as.numeric(paste(formula_split[(element_ind[i]+1):length(formula_split)], collapse = ""))
  # brokenFormula <- brokenFormula[which(brokenFormula$ElementNumber!=0),]
  return(brokenFormula[1,2])
}

get.c.diff <- function(index, network)
{
  c1 <- get.carbon.number(network$formula.1[index])
  c2 <- get.carbon.number(network$formula.2[index])
  c.diff <- abs(c1-c2)
  if (c.diff <= 6){
    return(c.diff)
  }
  else {
    return("carbon difference too large")
  }
}


#
# DS: some new functions for plotting
# we can move them to somewhere else
#
get_mid_list <- function(peak.x, peak.y, e, mi_data){
  # only for a single experiment
  x <- get_avg_mid(mi_data, peak.x, e)
  y <- get_avg_mid(mi_data, peak.y, e)

  if (length(x) == length(y)){
    mid.list <- list(x, y)
    name.list <- list(peak.x, peak.y)
  }
  else if (length(x) > length(y)){
    mid.list <- list(x, y, solution(x, y), find_convolution(x, y))
    name.list <- list(paste0("Peak: ", peak.x), paste0("Peak: ", peak.y), "Peak: unknown", paste0("unknown + ", peak.x))
  }
  else {
    mid.list <- list(y, x, solution(y, x), find_convolution(y, x))
    name.list <- list(paste0("Peak: ", peak.y), paste0("Peak: ", peak.x), "Peak: unknown", paste0("unknown + ", peak.x))
  }

  return(list(mid.list, name.list))
}


plot_mids <- function(mid.list, name.list, cosine, include_M0 = TRUE)
{
  # prevent R CMD CHECK warnings
  MI <- MID <- NULL

  p <- list()

  for (i in 1:length(mid.list)){

    mid <- mid.list[[i]]

    if (include_M0 == FALSE) {
      df <- data.frame(MI = 1:(length(mid)-1), MID = mid[2:length(mid)])
      df$MI <- factor(df$MI, levels = unique(df$MI))
    }
    else {
      df <- data.frame(MI = 0:(length(mid)-1), MID = mid)
      df$MI <- factor(df$MI, levels = unique(df$MI))
    }

    p[[i]] <- ggplot(df, aes(x = MI, y = MID)) +
      geom_bar(stat = "identity", width = 0.1, color = "white", fill = "skyblue") +
      theme_classic() + ylim(0, 1) +
      labs(x = "Mass Isotopomer (MI)", y = "MI fraction", title = name.list[[i]])
  }
  plot <- grid.arrange(grobs = p, ncol = 2, top = paste0("Cosine: ", cosine))
  return(plot)
}

