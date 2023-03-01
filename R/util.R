#
#  Some utility functions
#

#' remove any NA values in a vector x
remove_na <- function(x)
{
  return(x[!is.na(x)])
}


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
