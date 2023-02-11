#
#  Some utility functions
#

remove_na <- function(x)
{
  return(x[!is.na(x)])
}

# take maximum of x after discard missing values
# and return NA if there are no non-infinite values
# NOTE: I'm not convinced we should return NA in
# this case; see ?max for reasons for using +/-Inf
# Keeping this for compatibility for now
max_nonempty <- function(x)
{
  max_inf <- suppressWarnings(max(x, na.rm = TRUE))
  # replace -Inf by NA
  return(ifelse(is.infinite(max_inf), as.double(NA), max_inf))
}

min_nonempty <- function(x)
{
  min_inf <- suppressWarnings(min(x, na.rm = TRUE))
  # replace +Inf by NA
  return(ifelse(is.infinite(min_inf), as.double(NA), min_inf))
}
