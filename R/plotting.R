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
plot_matrix <- function(mat, limits = NULL)
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
      ggplot2::scale_fill_gradient(
        low = "white", high = "red",
        limits = limits,
        oob = scales::squish) +
      ggplot2::theme_bw()
  )
}


#' Plot a matrix of MIDs
#'
#' Here we discard MI 0 as it is redundant with the others.
#' To get axes labels right, the matrix should have dimnames set
#'
#' @param mid_mat A matrix whose columns are MID vectors
#' @param plot_title Title for the plot
#' @param max_mi_fraction MI fraction at which the heatmap saturates
#' @returns A ggplot object
#' @export
#'
plot_mid_matrix <- function(mid_mat, max_mi_fraction = 1.0, plot_title = "MID")
{
  return(
    plot_matrix(mid_mat[-1, ], limits = c(0, max_mi_fraction)) +
      labs(x = "Experiment", y = "MI", title = plot_title)
  )
}


#' Plot MID matrices for a convolution x * z = y
#'
#' @param mi_data An MIData object
#' @param index_x index of peak x
#' @param index_y index of peak y
#' @param index_z index of peak z
#' @returns A ggplot object
#' @export
#'
plot_convolution <- function(mi_data, index_x, index_y, index_z)
{
  # average MID matrices (MIs x experiments)
  mids_x <- get_avg_mid(mi_data, index_x)
  mids_y <- get_avg_mid(mi_data, index_y)
  mids_z <- get_avg_mid(mi_data, index_z)
  # check dimensions
  stopifnot(nrow(mids_x) + nrow(mids_z) - 2 == nrow(mids_y) - 1)
  # convolute x * y
  mids_xz <- sapply(
    1:length(mi_data$experiments),
    function(i) convolute(mids_x[, i], mids_z[, i])
  )
  # peak IDs
  peak_id_x <- mi_data$peak_ids[index_x]
  peak_id_y <- mi_data$peak_ids[index_y]
  peak_id_z <- mi_data$peak_ids[index_z]
  # construct plot and return it
  return(
    grid.arrange(
      plot_mid_matrix(mids_x, peak_id_x),
      plot_mid_matrix(mids_z, peak_id_z),
      plot_mid_matrix(mids_xz, paste(peak_id_x, peak_id_z, sep="*")),
      plot_mid_matrix(mids_y, peak_id_y),
      ncol = 2)
  )
}

