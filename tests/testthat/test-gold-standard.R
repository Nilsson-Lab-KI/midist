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


test_that("continuous_accuracy is correct", {
  # lower-diagonal elements
  dm_lower <- c(0.1, 0.5, 0.7, 0.3, 0.4, 0.9, 0.12, NA, 0.2, 0.2)
  gs_lower <- c(1, 0, 1, 1, 0, 0, 1, 0, 1, 0)
  # ordered elements
  dm_order <- c(1, 7, 9, 10, 4, 5, 2, 3, 6, 8)
  gs_ordered <- gs_lower[dm_order]

  accuracy_df <- continuous_accuracy(dm, gs)
  # recall = (no. true positives) / (no. true)
  expect_equal(
    accuracy_df$recall,
    cumsum(gs_ordered) / sum(gs_lower)
  )
  # precision = (no. true positives) / (no. positives)
  expect_equal(
    accuracy_df$precision,
    cumsum(gs_ordered) / 1:10
  )
})


test_that("get_continuous_accuracy is correct", {
  # number of true pairs
  gs_true = sum(gs)
  # true positive indicator in order of increasing distance
  # for off-diagonal elements of dm, with NA last and with ties in original order
  # dm_order <- c(1, 5, 4, 17, 12, 16, 19, 20, 3, 13, 7, 14, 2, 9, 6, 10, 8, 11, 15, 18)
  tp_ordered <- c(1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0)
  # list of 20 values, for all off-diagonal elements
  accuracy_df <- get_continuous_accuracy(
    dm, gs,
    # none of the remaining arguments are actually used by this function,
    # but are added to the returned data frame; this should be refactored
    measure = 1:20, subset_size = 1:20, subset_sample_no = 1:20,
    noise = 1:20, experiment = 1:20
  )
  # true positive rate = (no. true positives) / (no. true)
  expect_equal(
    accuracy_df$tpr,
    cumsum(tp_ordered) / gs_true
  )
  # precision = (no. true positives) / (no. positives)
  expect_equal(
    accuracy_df$precision,
    cumsum(tp_ordered) / 1:20
  )

})


test_that("get_global_percentile_accuracy is correct", {
  # percentiles of the 20 off-diagonal elements, ordered
  percentiles <- c(0, 10, 50, 80, 100)
  perc_index <- as.integer(quantile(1:20, percentiles / 100))

  accuracy <- get_global_percentile_accuracy(
    pairwise_matrix = dm,
    gold_standard = gs,
    percentiles = percentiles,
    # the remaining arguments are not used
    measure = 1:5, noise = 1:5, experiment = 1:5, subset_size = 1:5
  )
  # true positive rate = (no. true positives) / (no. true)
  expect_equal(
    accuracy$tpr,
    c(0.2, 0.2, 0.8, 0.8, 1.0)
  )
  # false positive rate = (no. false positives) / (no. false)
  expect_equal(
    accuracy$fpr,
    c(0.0, 0.0, 0.1333, 0.40, 0.5333),
    tolerance = 1e-3

  )

  # precision = (no. true positives) / (no. true)

  # false discovery rate = (no. false positives) / (no. positive)

})

