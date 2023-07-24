#
# Tests for function in distance.R
#


test_that("apply_no_m0 is correct", {
  x <- rnorm(n = 4)
  y <- rnorm(n = 4)
  f <- function(x, y) rlang::hash(c(x, y))
  expect_equal(
    apply_no_m0(f, x, y),
    f(x[-1], y[-1])
  )
})


test_that("cosine_sim is correct", {
  # identical vectors have similarity 1
  x <- rnorm(n = 4)
  expect_equal(cosine_sim(x, x), 1)
  # orthogonal vectors have similarity 0
  expect_equal(cosine_sim(c(1, 0), c(0, 1)), 0)
  # with a zero vector similarity is NA
  x <- rnorm(n = 4)
  y <- rep(0, 4)
  expect_true(is.na(cosine_sim(x, y)))
})


test_that("cosine_dist is correct", {
  # identical vectors have distance 0
  x <- rnorm(n = 4)
  expect_equal(cosine_dist(x, x), 0)
  # orthogonal vectors have distance 1
  expect_equal(cosine_dist(c(1, 0), c(0, 1)), 1)
  # with a zero vector distance is NA
  x <- rnorm(n = 4)
  y <- rep(0, 4)
  expect_true(is.na(cosine_dist(x, y)))
  # for any non-zero vectors, cosine_dist = 1 - cosine_sim
  x <- rnorm(n = 4)
  y <- rnorm(n = 4)
  expect_equal(cosine_dist(x, y), 1 - cosine_sim(x, y))
})


test_that("euclidean distance is correct", {
  # identical vectors have distance 0
  x <- rnorm(n = 4)
  expect_equal(euclidean_dist(x, x), 0)
  # distance between 2D unit vectors
  expect_equal(euclidean_dist(c(1, 0), c(0, 1)), sqrt(2))
  # with a zero vector, distance is equal to norm
  x <- rnorm(n = 4)
  y <- rep(0, 4)
  expect_equal(euclidean_dist(x, y), norm(x, type="2"))
})


test_that("squared euclidean distance is correct", {
  # identical vectors have similarity 0
  x <- rnorm(n = 4)
  expect_equal(euclidean_dist_sq(x, x), 0)
  # distance between 2D unit vectors
  expect_equal(euclidean_dist_sq(c(1, 0), c(0, 1)), 2)
  # with a zero vector, distance is equal to sum of squares
  x <- rnorm(n = 4)
  y <- rep(0, 4)
  expect_equal(euclidean_dist_sq(x, y), sum(x*x))
})


test_that("dot product similarity is correct", {
  # for identical vectors, similarity is sum of squares
  x <- rnorm(n = 4)
  expect_equal(dot_sim(x, x), sum(x*x))
  # similarity between unit vectors is zero
  expect_equal(dot_sim(c(1, 0), c(0, 1)), 0)
  # with a zero vector, similarity is zero
  x <- rnorm(n = 4)
  y <- rep(0, 4)
  expect_equal(dot_sim(x, y), 0)
})


test_that("dot product distance is correct", {
  # for identical vectors, distance is 0
  x <- rnorm(n = 4)
  expect_equal(dot_dist(x, x), 0)
  # similarity between unit vectors is 1
  expect_equal(dot_dist(c(1, 0), c(0, 1)), 1)
  # with a zero vector, distance is zero
  x <- rnorm(n = 4)
  y <- rep(0, 4)
  expect_equal(dot_dist(x, y), 0)
})


random_prob_vector <- function(n)
{
  x <- runif(n = 4)
  return(x / sum(x))
}


test_that("Jensen-Shannon distance is correct", {
  # two random probability vectors
  x <- random_prob_vector(n = 4)
  y <- random_prob_vector(n = 4)
  # distance is < 1 with probability 1
  expect_true(jensen_shannon(x, y) < 1)
  # if any vector component is zero, the result is NA
  x_zero <- x
  x_zero[1] <- 0
  x_zero <- x_zero / sum(x_zero)
  expect_true(is.na(jensen_shannon(x_zero, y)))
  # for identical vectors, distance is 0
  expect_equal(jensen_shannon(x, x), 0)
})


test_that("Manhattan distance is correct", {
  # identical vectors have distance 0
  x <- rnorm(n = 4)
  expect_equal(manhattan_distance(x, x), 0)
  # distance between unit vectors is 2
  expect_equal(manhattan_distance(c(1, 0), c(0, 1)), 2)
  # with a zero vector, distance is equal to sum
  x <- rnorm(n = 4)
  y <- rep(0, 4)
  expect_equal(manhattan_distance(x, y), sum(abs(x)))
})


test_that("Canberra distance is correct", {
  # identical vectors have distance 0
  x <- rnorm(n = 4)
  expect_equal(canberra_distance(x, x), 0)
  # distance between 2D unit vectors is 1/1 + 1/1 = 2
  expect_equal(canberra_distance(c(1, 0), c(0, 1)), 2)
  # with a zero vector, distance is equal to n
  x <- rnorm(n = 4)
  y <- rep(0, 4)
  expect_equal(canberra_distance(x, y), 4)
  # if any element is zero in both x and y the distance is NA
  x <- rnorm(n = 4)
  x[1] <- 0
  y <- rnorm(n = 4)
  y[1] <- 0
  expect_true(is.na(canberra_distance(x, y)))
})


test_that("Bray-Curtis distance is correct", {
  # identical vectors have distance 0
  x <- rnorm(n = 4)
  expect_equal(bray_curtis_distance(x, x), 0)
  # distance between unit vectors is 1
  expect_equal(bray_curtis_distance(c(1, 0), c(0, 1)), 1)
  # with one zero vector, distance is 1
  x <- rnorm(n = 4)
  y <- rep(0, 4)
  expect_equal(bray_curtis_distance(x, y), 1)
  # with two zero vectors, distance is NA
  expect_true(is.na(bray_curtis_distance(y, y)))
})


test_that("mi_weighted_distance is correct", {
  # identical vectors have distance 0
  x <- rnorm(n = 4)
  expect_equal(mi_weighted_distance(x, x), 0)
  # distance between unit vectors depends on dimension
  for(n in 1:5) {
    x <- c(1, rep(0, n))
    y <- c(rep(0, n), 1)
    expect_equal(mi_weighted_distance(x, y), n)
  }
})


test_that("mi_weighted_dist_normalized is correct", {
  # identical vectors have distance 0
  x <- rnorm(n = 4)
  expect_equal(mi_weighted_dist_normalized(x, x), 0)
  # distance between unit vectors is always 1
  for(n in 1:5) {
    x <- c(1, rep(0, n))
    y <- c(rep(0, n), 1)
    expect_equal(mi_weighted_dist_normalized(x, y), 1)
  }
})


