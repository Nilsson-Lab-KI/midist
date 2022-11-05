
test_that("cosine_sim is correct", {
  # identical vectors have similarity 1
  expect_equal(cosine_sim(c(2,3), c(2,3)), 1)
  # orthogonal vectors have zero similarity
  expect_equal(cosine_sim(c(1,0), c(0,1)), 0)
})
