#
#  Some utility functions
#


#' Take maximum of x after discarding missing values
#'
#' @param x A vector of values. If x has an "index"
#' attribute, it must be a vector of same length as x.
#'
#' @returns The largest value in x, or NA if x is empty. If x has an "index"
#' attribute, the element corresponding to the largest value is added as
#' an "index" attribute to the result.
#' @export
max_nonempty <- function(x)
{
  best_nonempty(x, function(x) order(x, decreasing = TRUE)[1])
}


#' Take minimum of x after discarding missing values
#'
#' @param x A vector of values. If x has an "index"
#' attribute, it must be a vector of same length as x.
#'
#' @returns The smallest value in x, or NA if x is empty. If x has an "index"
#' attribute, the element corresponding to the smallest value is added as
#' an "index" attribute to the result.
#' @export
min_nonempty <- function(x)
{
  best_nonempty(x, function(x) order(x, decreasing = FALSE)[1])
}


#' Apply a selection function to x while handling the "index" attribute
#'
#' NOTE: I'm not convinced we should return NA in
#' case of empty lists; see ?max for reasons for using +/-Inf
#' Keeping this for compatibility for now
#'
#' @param x A vector of values. If x has an "index"
#' attribute, it must be a vector of same length as x.
#' @param select A function selecting a value from x, returning its index
#'
#' @returns The value in x chosen by select, or NA if x is empty. If x has an "index"
#' attribute, the element corresponding to the selected value is added as
#' an "index" attribute to the result.
#' @export
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


#' Get object without attributes
#' @param x Any object
#' @returns The object x with attributes(x) == NULL
without_attr <- function(x)
{
  x_no_attr <- x
  attributes(x_no_attr) <- NULL
  return(x_no_attr)
}


#' Get object with an attribute added
#' @param x Any object
#' @param attr_name name of the attribut to add
#' @param attr_value the attribute value to add
#' @returns The object x with attr(x, attr_name) == attr_value
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


# utility for handling multiple assignments
# taken from
# (originally called list[...] but I felt that gets confused with base::list)

#' Unique object used to assign to list elements
#' @export
assign_list <- structure(NA, class = "AssignList")


#' Operator for assignment to an AssignList
#'
#' For example, assign_list[x, y] <- list(a,b) evaluates to
#' "[<-AssignList"(assign_list, x, y, some_list) and is effectively the same as
#' x <- a and y <- b.
#' See help for "[<-"]
#
"[<-.AssignList" <- function(assign_list, ..., value)
{
  # get arguments (first element is function name)
  args <- as.list(match.call())
  # take arguments after assign_list, before value = elements to assign
  args <- args[-c(1:2, length(args))]
  #
  length(value) <- length(args)
  for(i in seq(along = args)) {
    arg <- args[[i]]
    if(!missing(arg))
      eval.parent(substitute(
        a <- v, list(a = arg, v = value[[i]])))
  }
  return(assign_list)
}


#' Test if a matrix is a valid distance matrix (metric)
#' @param mat A (candidate) distance matrix
#' NOTE: this function should be fixed, see https://github.com/Nilsson-Lab-KI/remn/issues/40
is_distance_matrix <- function(mat) {
  # necessary intervention to prevent the floating number error
  mat <- as.matrix(Matrix::nearPD(mat)$mat)
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
  }

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
