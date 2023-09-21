#
# Some plotting functions for mass isotopomer distributions
#


#' Plot MID as barchart with error bars
#'
#' @param mid_mat A matrix whose columns are MID vectors
#' @returns A ggplot object
#' @importFrom stats sd
#' @export
#'
plot_mid_barchart <- function(mid_mat)
{
  n <- nrow(mid_mat) - 1
  # compute mean and sd for each MI fraction
  mid_df <- data.frame(
    mi = 0:n,
    mean = rowMeans(mid_mat),
    sd = apply(mid_mat, 1, stats::sd)
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


#' Plot a matrix as a heatmap
#'
#' @param mat Any real matrix

#' @importFrom reshape2 melt
#' @returns A ggplot object
#'
plot_matrix <- function(mat)
{
  mat_melted <- reshape2::melt(mat)
  return(
    ggplot2::ggplot(
      mat_melted,
      ggplot2::aes(
        x = mat_melted[, 2],
        y = mat_melted[, 1])
    ) +
      ggplot2::geom_raster(aes(fill = value)) +
      ggplot2::scale_fill_gradient(low = "white", high = "red") +
      ggplot2::theme_bw()
  )
}


#' Plot a matrix of MIDs
#'
#' To get axes labels right, the matrix should have dimnames set
#'
#' @param mid_mat A matrix whose columns are MID vectors
#' @param plot_title Title for the plot
#' @returns A ggplot object
#' @export
#'
plot_mid_matrix <- function(mid_mat, plot_title = "MID")
{
  return(
    plot_matrix(mid_mat) +
      labs(x = "Experiment", y = "MI", title = plot_title)
  )
}


