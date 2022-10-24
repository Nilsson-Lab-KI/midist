
#
# MID normalization
#

test_that("normalized MIDs sum to 1", {
  # example "peak area" data, each column an MID
  peak_areas <- matrix(
    c(1,0,
      3,2,
      7,0),
    nrow = 3, ncol = 2)
  mids <- normalize_mids(peak_areas)
  expect_equal(colSums(mids), c(1,1))
})

test_that("normalizing zero vector gives zero vector", {
  # example "peak area" data, each column an MID
  peak_areas <- matrix(
    c(0,
      0,
      0),
    nrow = 3, ncol = 1)
  mids <- normalize_mids(peak_areas)
  expect_equal(colSums(mids), c(0))
})
