#
# Tests for "utility" function in util.R
#

test_that("max_nonempty is correct", {
  # full vector without index
  a <- c(0.1, 0.3, 0.7, 0.0)
  max_a <- max_nonempty(a)
  expect_equal(max_a, with_attr(0.7, "index", as.integer(NA)))
  expect_true(max_a == 0.7)

  # vector with some NAs
  a <- c(0.1, NA, 0.7, NA)
  max_a <- max_nonempty(a)
  expect_equal(max_a, with_attr(0.7, "index", as.integer(NA)))

  # vector with all NAs
  a <- c(NA, NA, NA)
  max_a <- max_nonempty(a)
  expect_true(is.na(max_a))

  # vector with index attribute
  a <- with_attr(c(0.1, NA, 0.7, NA), "index", c(1, NA, 2, NA))
  max_a <- max_nonempty(a)
  # expected result with attribute
  expected <- 0.7
  attr(expected, "index") <- 2
  expect_equal(max_a, expected)
  expect_true(max_a == 0.7)

  # empty vector yields NA
  expect_true(is.na(max_nonempty(c())))
})


test_that("min_nonempty is correct", {
  # full vector without index
  a <- c(0.1, 0.3, 0.7, 0.05)
  min_a <- min_nonempty(a)
  expect_equal(min_a, with_attr(0.05, "index", as.integer(NA)))
  expect_true(min_a == 0.05)

  # vector with some NAs
  a <- c(0.1, NA, 0.7, NA)
  min_a <- min_nonempty(a)
  expect_equal(min_a, with_attr(0.1, "index", as.integer(NA)))

  # vector with all NAs
  a <- c(NA, NA, NA)
  min_a <- min_nonempty(a)
  expect_true(is.na(min_a))

  # vector with index attribute
  a <- with_attr(c(0.1, NA, 0.7, NA), "index", c(1, NA, 2, NA))
  min_a <- min_nonempty(a)
  # expected result with attribute
  expected <- 0.1
  attr(expected, "index") <- 1
  expect_equal(min_a, expected)
  expect_true(min_a == 0.1)

  # empty vector yields NA
  expect_true(is.na(min_nonempty(c())))
})


