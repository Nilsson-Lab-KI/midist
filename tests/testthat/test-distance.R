
test_that("cosine_sim is correct", {
  # identical vectors have similarity 1
  expect_equal(cosine_sim(c(2,3), c(2,3)), 1)
  # orthogonal vectors have zero similarity
  expect_equal(cosine_sim(c(1,0), c(0,1)), 0)
  # with a zero vector similarity is NA
  expect_true(is.na(cosine_sim(c(1,0), c(0,0))))
})

test_that("cosine_dist is correct", {
  # identical vectors have similarity 0
  expect_equal(cosine_dist(c(2,3), c(2,3)), 0)
  # orthogonal vectors have distance 1
  expect_equal(cosine_dist(c(1,0), c(0,1)), 1)
  # with a zero vector distance is NA
  expect_true(is.na(cosine_dist(c(1,0), c(0,0))))
})
