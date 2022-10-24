#
# test MID handling functions
#

test_that("enrichment values are correct", {
  # unit vector has enrichment zero
  expect_equal(isotopic_enrichment(c(1,0,0)), 0)
  # fully labeled vector has enrichment 1
  expect_equal(isotopic_enrichment(c(0,0,1)), 1)
  # 50/50 vector has enrichment 0.5
  expect_equal(isotopic_enrichment(c(0.5,0,0.5)), 0.5)
})

test_that("convolution matrix has right dimensions", {
  # convolution 2-carbon x with 3-carbon y
  A <- convolution_matrix(x = c(0.1,0.3,0.6), y_carbons = 3)
  expect_equal(dim(A), c(2+3+1, 3+1))
})

test_that("MID convolution is correct", {
  # two 1-carbon unit vectors
  expect_equal(convolute(c(1,0), c(1,0)), c(1,0,0))
  # short * long vector
  expect_equal(convolute(c(0.5, 0.5), c(0.2, 0, 0.8)), c(0.1, 0.1, 0.4, 0.4))
  # long * short vector
  expect_equal(convolute(c(0.2, 0, 0.8), c(0.5, 0.5)), c(0.1, 0.1, 0.4, 0.4))

})
