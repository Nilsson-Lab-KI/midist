#
# Unit tests for plottning functions.
#

mid_mat <- matrix(
  c(
    0.10, 0.15,
    0.05, 0.10,
    0.80, 0.70,
    0.05, 0.05
  ),
  nrow = 4,
  dimnames = list(mi = 0:3, experiments = 1:2)
)

test_that("plot_mid_barchart yields a ggplot object", {
  expect_s3_class(
    plot_matrix(mid_mat),
    "ggplot"
  )
})

test_that("plot_mid_matrix yields a ggplot object", {
  expect_s3_class(
    plot_mid_matrix(mid_mat),
    "ggplot"
  )
})


