
test_that("get_continuous_accuracy is correct", {
  # example 5 x 5 distance matrix
  dm <- matrix(
    c(
      0,    0.1,  0.5,  0.3,  0.12,
      0.1,  0,    0.7,  0.4,  NA,
      0.5,  0.7,  0,    0.9,  0.2,
      0.3,  0.4,  0.9,  0,    0.2,
      0.12, NA, 0.2, 0.2,   0
    ),
    nrow = 5
  )
  # gold standard (true pair = 1)
  gs <- matrix(
    c(
      0, 1, 0, 1, 1,
      1, 0, 1, 0, 0,
      0, 1, 0, 0, 1,
      1, 0, 0, 0, 0,
      1, 0, 1, 0, 0
    ),
    nrow = 5
  )
  gs_positive = sum(gs)
  # true positive indicator in order of increasing distance
  # for off-diagonal elements of dm, with NA last and with ties in original order
  # dm_order <- c(1, 5, 4, 17, 12, 16, 19, 20, 3, 13, 7, 14, 2, 9, 6, 10, 8, 11, 15, 18)
  tp_ordered <- c(1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0)
  # list of 20 values, for all off-diagonal elements
  accucary_df <- get_continuous_accuracy(
    dm, gs,
    # none of the remaining arguments are actually used by this function,
    # but are added to the returned data frame; this should be refactored
    measure = 1:20, subset_size = 1:20, subset_sample_no = 1:20,
    noise = 1:20, noise_sample_no = 1:20, experiment = 1:20
  )
  # true positive rate = (no. true positives) / (no. true)
  expect_equal(
    accucary_df$tpr,
    cumsum(tp_ordered) / gs_positive
  )
  # precision = (no. true positives) / (no. true)
  expect_equal(
    accucary_df$precision,
    cumsum(tp_ordered) / 1:20
  )

})
