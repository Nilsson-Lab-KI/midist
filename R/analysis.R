# 
# functions to analyse similarity matrices
#

# reads similarity matrix and fixes row/column names
read_similarity_matrix <- function(file_name, file_path)
{
  sm <- as.data.frame(read.delim(file.path(file_path, file_name), 
                                 header = T, sep = "\t"))
  
  # reading from file, column names get messed up, so they need a fix
  rownames(sm) <- colnames(sm) <- sub("\\.", "-", sub("X", "", colnames(sm)))
  
  return(sm)
}

# fixes missing values by assigning ones or zeros to them depending on whether minimal is T of F, respectively
fix_na <- function(similarity_matrix, minimal = T)
{
  similarity_matrix[which(is.na(similarity_matrix), arr.ind = T)] <- 0
  if (minimal == F)
    similarity_matrix[which(is.na(similarity_matrix), arr.ind = T)] <- 1
  
  return(similarity_matrix)
}

# make cluster trees for the full sets
cluster_tree <- function(similarity_matrix, distance = F)
{
  
  distance_matrix <- 1 - similarity_matrix
  
  if (distance == T)
    distance_matrix <- similarity_matrix
  
  # the clustering tree
  dendro_labels <- colnames(distance_matrix)
  
  return(hclust(as.dist(distance_matrix)))
}



