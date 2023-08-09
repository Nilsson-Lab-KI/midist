#
# test MID handling functions
#

test_that("isotopic_enrichment is correct", {
  # unit vector has enrichment zero
  expect_equal(isotopic_enrichment(c(1, 0, 0)), 0)
  # fully labeled vector has enrichment 1
  expect_equal(isotopic_enrichment(c(0, 0, 1)), 1)
  # 50/50 vector has enrichment 0.5
  expect_equal(isotopic_enrichment(c(0.5, 0, 0.5)), 0.5)
})

test_that("filter_enrichment works correctly", {
  expect_equal(
    filter_enrichment(c(0.99, 0.01), tol = 0.01),
    c(0, 0)
  )
  expect_equal(
    filter_enrichment(c(0.99, 0.01), tol = 0.01 - 1e-10),
    c(0.99, 0.01)
  )
})

test_that("add_gamma_noise works correctly", {
  # the noisy MID should sum to 1
  expect_equal(sum(add_gamma_noise(c(0.8, 0.01, 0.1, 0.09), stdev = 0.01)), 1)
  # a zero vector results in NA
  expect_equal(is.na(max(add_gamma_noise(c(0,0,0,0), stdev = 0.01))), T)
  # unlabelled remains unlabelled
  expect_equal(add_gamma_noise(c(1,0,0,0), stdev = 5), c(1,0,0,0))
})

test_that("convolution_matrix gives correct dimensions", {
  # convolution 2-carbon x with 3-carbon y
  A <- convolution_matrix(x = c(0.1, 0.3, 0.6), y_carbons = 3)
  expect_equal(dim(A), c(2 + 3 + 1, 3 + 1))
})

test_that("convolute is correct", {
  # two 1-carbon unit vectors
  expect_equal(convolute(c(1, 0), c(1, 0)), c(1, 0, 0))
  # short * long vector
  expect_equal(convolute(c(0.5, 0.5), c(0.2, 0, 0.8)), c(0.1, 0.1, 0.4, 0.4))
  # long * short vector
  expect_equal(convolute(c(0.2, 0, 0.8), c(0.5, 0.5)), c(0.1, 0.1, 0.4, 0.4))
  # convolutions with zero vector
  expect_equal(convolute(c(1, 0), c(0, 0)), c(0, 0, 0))
  expect_equal(convolute(c(0, 0), c(1, 0)), c(0, 0, 0))
})

test_that("convolute_cols agrees with convolute", {
  x <- c(0.8, 0.2)
  # random MID matrix, each column is one MID
  y_mat <- matrix(runif(3 * 4), nrow = 3, ncol = 4)
  y_mat <- t(t(y_mat) / colSums(y_mat))
  # test convolution
  expect_equal(
    convolute_cols(x, y_mat),
    apply(y_mat, MARGIN = 2, function(y) convolute(x, y))
  )
})


test_that("c13correct is correct", {
  # case where simplex constraints must be enforced
  expect_equal(c13correct(c(1, 0, 0)), c(1, 0, 0))
  # case from simplex interior to corner
  expect_equal(c13correct(c(0.9, 0.1), p = 0.1), c(1, 0))
  # case within interior of simplex
  expect_equal(c13correct(c(0.1, 0.9), p = 0.2), c(0.125, 0.875))
  # edge case for zero atoms
  expect_equal(c13correct(c(1)), c(1))
})


test_that("solution is correct", {
  # find y minimizing  || x * y - z ||
  x <- c(0.5, 0.5)
  y <- c(0.2, 0, 0.8)
  z <- c(0.1, 0.1, 0.4, 0.4)
  expect_equal(solution(z, x), y)
})


test_that("find_convolution is correct", {
  # find y minimizing  || x * y - z ||
  x <- c(0.5, 0.5)
  #y <- c(0.2, 0, 0.8)
  z <- c(0.1, 0.1, 0.4, 0.4)
  expect_equal(find_convolution(z, x), z)
})


test_that("random_mid is correct", {
  mean_mid <- c(0.1, 0.2, 0.0, 0.7)
  n <- 10
  random_mids <- random_mid(mean_mid, 0.1, n)
  # check dimensions
  expect_true(is.matrix(random_mids))
  expect_equal(
    dim(random_mids),
    c(4, n)
  )
  # all MI fractions must be positive
  expect_true(all(random_mids >= 0))
  # zeros in the mean vector gives zeros in all observations
  expect_true(all(random_mids[3, ] == 0))
  # MIDs must sum to 1
  expect_equal(
    colSums(random_mids),
    rep(1.0, n)
  )
  # with n == 1 we obtain a matrix with 1 column
  random_mids <- random_mid(mean_mid, 0.1, 1)
  expect_true(is.matrix(random_mids))
  expect_equal(dim(random_mids), c(4, 1))
})

