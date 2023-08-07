#
#  Some utility functions
#


#' Take maximum of x after discarding missing values,
#' return NA if the result is empty. If x has an "index"
#' attribute, it must be a list of same length as x;
#' the element corresponding to the maximum is added as
#' an "index" attribute to the result
#'
#' NOTE: I'm not convinced we should return NA in
#' case of empty lists; see ?max for reasons for using +/-Inf
#' Keeping this for compatibility for now
max_nonempty <- function(x)
{
  best_nonempty(x, function(x) order(x, decreasing = TRUE)[1])
}


#' Similar to max_nonempty, but takes the minimum
#' @export
min_nonempty <- function(x)
{
  best_nonempty(x, function(x) order(x, decreasing = FALSE)[1])
}


# generic function taking a "select" function picking the "best" element
best_nonempty <- function(x, select)
{
  if(length(x) == 0)
    return(with_attr(as.double(NA), "index", as.integer(NA)))
  # find maximum index and value (possibly NA)
  x_index <- select(x)
  best_value <- x[x_index]
  # handle index attribute
  best_index <- ifelse(is.null(attr(x, "index")),
                       as.integer(NA),
                       attr(x, "index")[x_index])
  return(with_attr(best_value, "index", best_index))
}


#' get object without attributes
without_attr <- function(x)
{
  x_no_attr <- x
  attributes(x_no_attr) <- NULL
  return(x_no_attr)
}


#' get object with an attribute added
with_attr <- function(x, attr_name, attr_value)
{
  x_attr <- x
  attr(x_attr, attr_name) <- attr_value
  return(x_attr)
}


#' Apply a function on MIDs x,y after removing M+0
#'
#' This simply removes the first component of each vector before applying f
#' @param f function to apply
#' @param x a vector
#' @param y a vector
#' @returns the result of calling f
#' @export
#'
apply_no_m0 <- function(f, x, y) {
  return(f(x[-1], y[-1]))
}


#' Test if a matrix is a valid distance matrix (metric)
#' NOTE: this function should be fixed, see https://github.com/Nilsson-Lab-KI/remn/issues/40
is_distance_matrix <- function(mat) {
  # necessary intervention to prevent the floating number error
  mat <- as.matrix(nearPD(mat)$mat)
  # Initialize a character vector to store failure reasons
  failures <- character(0)

  # Check if the input is a matrix
  if (!is.matrix(mat)) {
    failures <- append(failures, "Input is not a matrix")
  }

  # Check if the matrix is square
  if (nrow(mat) != ncol(mat)) {
    failures <- append(failures, "Matrix is not square")
  }

  # Check if the matrix is symmetric
  if (!isSymmetric(mat)) {
    failures <- append(failures, "Matrix is not symmetric")
  }

  # Check if the diagonal entries are zero
  if (!all(diag(mat) == 0)) {
    failures <- append(failures, "Diagonal entries are not zero")
  }

  # Check if the matrix is non-negative
  if (!all(mat >= 0)) {
    failures <- append(failures, "Matrix contains negative values")
  }R

  # Check the triangle inequality
  for (i in 1:nrow(mat)) {
    for (j in 1:nrow(mat)) {
      for (k in 1:nrow(mat)) {
        if (round(mat[i,j], digits = 3) + round(mat[j,k], digits = 3) < round(mat[i,k], digits = 3)) {
          failures <- append(failures, "Triangle inequality fails")
        }
      }
    }
  }

  # If no failures, return NULL, else return the failure reasons
  if (length(failures) == 0) {
    return(NULL)
  } else {
    return(failures)
  }
}
