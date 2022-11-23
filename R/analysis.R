# 
# functions to analyse similarity matrices
#

# reads similarity matrix and fixes row/column names
#' @export
read_similarity_matrix <- function(file_name, file_path)
{
  sm <- as.data.frame(read.delim(file.path(file_path, file_name), 
                                 header = T, sep = "\t"))
  
  # reading from file, column names get messed up, so they need a fix
  rownames(sm) <- colnames(sm) <- sub("\\.", "-", sub("X", "", colnames(sm)))
  
  return(sm)
}

# fixes missing values by assigning ones or zeros to them depending on whether minimal is T of F, respectively
#' @export
fix_na <- function(similarity_matrix, minimal = T)
{
  similarity_matrix[which(is.na(similarity_matrix), arr.ind = T)] <- 0
  if (minimal == F)
    similarity_matrix[which(is.na(similarity_matrix), arr.ind = T)] <- 1
  
  return(similarity_matrix)
}

#' @export
distribution_plot <- function(similarity_matrix, color = "black")
{
  # convert similarity matrix into a sorted vector
  sm_v <- as.vector(similarity_matrix)[order(as.vector(similarity_matrix))]
  # create a data frame for plotting
  sm_df <- data.frame(pair = 1:length(sm_v),
                      similarity = sm_v)
  # plot the distribution
  distribution_plot <- ggplot(sm_df, aes(pair, similarity)) + geom_point(color = color) +
    theme_classic() + xlab("Metabolite pair") + ylab("Similarity")
  
  return(distribution_plot)
}

#' @export
clustered_heatmap <- function(similarity_matrix)
{
  # convert similarities into distances
  dm <- 1 - similarity_matrix
  
  # basic heatmap
  clustered_heatmap <- heatmap(dm, scale = "none")
  
  return(clustered_heatmap)
}

#' @export
interactive_heatmap <- function(similarity_matrix)
{
  # convert similarities into distances
  dm <- 1 - similarity_matrix
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
cluster_tree <- function(similarity_matrix, method = "average")
{
  # convert similarities into distances
  dm <- 1 - similarity_matrix
  # metabolite names / tree labels
  dendro_labels <- colnames(dm)
  # the clustering tree
  tree <- hclust(as.dist(dm), method = method)
  
  return(tree)
}

#' @export
scatter <- function(matrix_list, color = "black")
{
  # convert these matrices into rows of a data frame
  df <- data.frame(similarity_matrix = as.numeric(matrix_list[[1]]),
                   gold_standard = as.numeric(matrix_list[[2]]))
  scatter_plot <- ggplot(df, aes(similarity_matrix, gold_standard)) + geom_point(color = color) +
    theme_classic() + xlab("Similarity matrix") + ylab("Gold standard")
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


