#
# Some plotting functions for mass isotopomer distributions
#


#' Plot MID as barchart with error bars
#'
#' @param mid_mat A matrix whose columns are MID vectors
#' @returns A ggplot object
#' @export

plot_mid_barchart <- function(mid_mat)
{
  n <- nrow(mid_mat) - 1
  mid_df <- data.frame(
    mi = 0:n,
    mean = rowMeans(mid_mat),
    sd = apply(mid_mat, 1, sd)
  )
  return(
    ggplot2::ggplot(mid_df) +
      ggplot2::geom_bar(
        ggplot2::aes(x = mi, y = mean),
        stat = "identity"
      ) +
      ggplot2::geom_errorbar(
        ggplot2::aes(
          x = mi,
          ymin = mean - sd,
          ymax = mean + sd
        ),
        width = 0.5
      ) +
      ggplot2::theme_classic()
  )
}


