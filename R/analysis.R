# 
# functions to analyse similarity matrices
#

# reads similarity matrix and fixes row/column names
#' @export
read_matrix <- function(file_name, file_path)
{
  matri <- as.data.frame(read.delim(file.path(file_path, file_name), 
                                 header = T, sep = "\t", check.names = F))
  
  return(matri)
}

# fixes missing values by assigning a specified value
#' @export
fix_na <- function(matri, value)
{
  matri[which(is.na(matri), arr.ind = T)] <- value
  return(matri)
}


#' @export
scatter <- function(results_df, c13_correction = F, measure = "cosine", 
                    color = "black", marginals = "histogram", 
                    cor_method = "spearman")
{
  # subset results_df based on input
  df <- results_df[which(results_df$c13_correction == c13_correction & results_df$measure_name == measure),]
  # remove NAs and infs from data frame
  df <- df[-which(is.na(df$distance) | is.infinite(df$distance)),]
  
  # calculate correlation
  corr <- cor(df$distance, df$gs_fraction, method = cor_method)
  
  # plot name
  if (c13_correction == F)
    plot_title <- paste("Measure:", measure, "distance | 13C correction: False") else
      plot_title <- paste("Measure:", measure, "distance | 13C correction: True")
  
  # standard scatter :
  p <- ggplot(df, aes(x = distance, y = gs_fraction, text = text)) +
    # set color and add transparency for a better density visualization
    geom_point(color = color, alpha = 0.1) +
    # turn off legend and remove the grid
    theme(legend.position="none") + theme_classic() +
    # set plot axis titles
    ggtitle(plot_title,
            subtitle = paste(cor_method, "correlation =", format(corr, digits = 2))) +
    xlab("Distance") + ylab("Fraction in gold standard")
  
  # with marginal histogram
  if (is.null(marginals) == F)
    p <- ggMarginal(p, type = marginals)
  
  return(p)
  
}


#' @export
distribution_plot <- function(matri, color = "black", measure_type = "Similarity")
{
  # convert similarity matrix into a sorted vector
  v <- as.vector(matri)[order(as.vector(matri))]
  # create a data frame for plotting
  df <- data.frame(pair = 1:length(v),
                      measure = v)
  # plot the distribution
  distribution_plot <- ggplot(df, aes(pair, measure)) + geom_point(color = color) +
    theme_classic() + xlab("Metabolite pair") + ylab(measure_type)
  
  return(distribution_plot)
}

#' @export
clustered_heatmap <- function(matri, distance = F)
{
  # convert similarities into distances
  dm <- 1 - matri
  
  if (distance == T)
    dm <- matri
  
  # basic heatmap
  clustered_heatmap <- heatmap(dm, scale = "none")
  
  return(clustered_heatmap)
}

#' @export
interactive_heatmap <- function(matri, distance = F)
{
  # convert similarities into distances
  dm <- 1 - matri
  
  if (distance == T)
    dm <- matri
  
  # create a distance data frame for plotting
  df_dm <- cbind(expand.grid(rownames(dm), rownames(dm)), expand.grid(dm))
  colnames(df_dm) <- c("Metabolite_1", "Metabolite_2", "Distance")
  
  # plot an interactive heatmap
  # keep in mind that this is not clustered
  interactive_heatmap <- plotly::plot_ly(
    data = df_dm,
    x = ~Metabolite_1, y = ~Metabolite_2, z = ~Distance, text = ~paste0('Metabolites: ', Metabolite_1, " & ", Metabolite_2),
    hoverinfo = "text",
    type = "heatmap"
  )
  
  return(interactive_heatmap)
}

#' @export
cluster_tree <- function(matri, distance = F, method = "average")
{
  # convert similarities into distances
  dm <- 1 - matri
  
  if (distance == T)
    dm <- matri
  
  # metabolite names / tree labels
  dendro_labels <- colnames(dm)
  # the clustering tree
  tree <- hclust(as.dist(dm), method = method)
  
  return(tree)
}

#' @export
scatter_old <- function(matrix_list, color = "black", correlation_method = "spearman",
                    x_name = "Similarity matrix", y_name = "Gold standard")
{
  # convert these matrices into rows of a data frame
  df <- data.frame(measure = as.numeric(matrix_list[[1]]),
                   gold_standard = as.numeric(matrix_list[[2]]))
  scatter_plot <- ggplot(df, aes(measure, gold_standard)) + geom_point(color = color) +
    theme_classic() + xlab(x_name) + ylab(y_name) + 
    ggtitle(paste0("Spearman corr: ", 
                   format(cor(df$measure, df$gold_standard, method = correlation_method), digits = 2)))
  return(scatter_plot)
}



#' @export
subset_matrix <- function(carbon_group, matrix, midata)
{
  peaks <- midata$n_atoms_index[match(carbon_group, names(midata$n_atoms_index))][[1]]
  return(matrix[peaks, peaks])
}

#' @export
intersection_matrix <- function(matrix_1, matrix_2)
{
  
  # create a list to store both matrices
  output <- list()
  mets <- intersect(colnames(matrix_1), colnames(matrix_2))
  
  output[[1]] <- as.matrix(matrix_1[match(mets, colnames(matrix_1)),
                          match(mets, colnames(matrix_1))], drop = F)
  
  output[[2]] <- as.matrix(matrix_2[match(mets, colnames(matrix_2)),
                          match(mets, colnames(matrix_2))], drop = F)
  
  return(output)
}

#' @export
add_text <- function(row){
  return(paste("Met_1:",as.character(row[1]),"\n",
               "Met_2:",as.character(row[2]),"\n",
               "Measure:",as.character(row[3]),"distance","\n"))
}


#' @export
analyse <- function(pairwise_matrix, midata, gold_standard = NULL, distance = F){
  pm <- list()
  class(pm) <- "PairwiseMatrix"
  
  pm$pairwise_matrix <- pairwise_matrix
  
  # plot the overall similarity / distance distribution
  pm$distribution_plot <- distribution_plot(pm$pairwise_matrix)
  
  # plot a clustered heatmap 
  pm$clustered_heatmap <- clustered_heatmap(pm$pairwise_matrix, distance = distance)
  
  # plot hierarchical clustering tree
  pm$cluster_tree <- cluster_tree(pm$pairwise_matrix, distance = distance)
  
  # plot an interactive heatmap - NOT CLUSTERED
  pm$interactive_heatmap <- interactive_heatmap(pm$pairwise_matrix, distance = distance)
  
  # separate C analysis
  carbon_groups <- as.vector(as.data.frame(table(input_data$midata$peak_n_atoms))[-which(as.data.frame(table(input_data$midata$peak_n_atoms))$Freq == 1),1])
  carbon_subsets <- lapply(carbon_groups, 
                           subset_matrix, 
                           pm$pairwise_matrix,
                           input_data$midata)
  
  # distribution plots for carbon groups
  pm$carbon_distributions <- lapply(carbon_subsets, distribution_plot)
  # clustered heatmaps for carbon groups
  pm$clustered_carbon_heatmaps <- lapply(carbon_subsets, clustered_heatmap, distance = distance)
  # hierarchical clustering trees for carbon groups
  pm$carbon_cluster_tree <- lapply(carbon_subsets, cluster_tree, method = "average", distance = distance)
  # interactive non clustered heatmaps for carbon groups
  pm$interactive_carbon_heatmaps <- lapply(carbon_subsets, interactive_heatmap, distance = distance)
  
  # name the list elements so that they are consistent with the midata
  names(pm$carbon_distributions) <- 
    names(pm$clustered_carbon_heatmaps) <- 
    names(pm$carbon_cluster_tree) <- 
    names(pm$interactive_carbon_heatmaps) <- as.character(carbon_groups)
  
  
  # compare similarity matrix to gold standard if there is a gold_standard input
  if (is.null(gold_standard) == F){
    
    if (is.character(gold_standard)) # read gold standard from file
      gold_standard <- as.data.frame(read.delim(gold_standard, header = T, sep = "\t", check.names = F))
    stopifnot(is.data.frame(gold_standard) == T | is.matrix(gold_standard) == T)
    
    if (dim(gold_standard)[2] - dim(gold_standard)[1] == 1)
      gold_standard <- gold_standard[,-1]
    pm$gold_standard <- gold_standard
    
    # all vs all scatter plot 
    # keep in mind that scatter() is defined by remn and not imported from another package
    common_met_matrices <- intersection_matrix(pm$pairwise_matrix, pm$gold_standard)
    pm$compare_to_goldstandard <- scatter(common_met_matrices)
    
    # scatter plots for each carbon group
    carbon_common_met_matrices <- lapply(carbon_subsets, intersection_matrix, pm$gold_standard)
    pm$carbon_compare_to_goldstandard <- lapply(carbon_common_met_matrices, scatter)
    names(pm$carbon_compare_to_goldstandard) <- as.character(carbon_groups)
    
  }
  
  return(pm)
  
}
