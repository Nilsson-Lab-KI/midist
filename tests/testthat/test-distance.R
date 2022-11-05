
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


test_that("conv_similarity returns a valid similarity score", {
  # peak area data, one experiment
  peak_areas = data.frame(
    peak_id <- c(rep("a",2), rep("b",3), rep("c",3), rep("d",4)),
    exp1 <- c(
      0.98, 0.02,
      0.8, 0.1, 0.1,
      0.6, 0.2, 0.2,
      0.8, 0.05, 0.1, 0.05))
  # create MIData object
  midata <- MIData(peak_areas, exp_names = "exp1")
  # for a vs d, we take maximum of a*b vs d and a*c vs d
  expect_equal(
    conv_similarity(midata, 1, 4, 1, cosine_sim), 0.9949405, tolerance = 1e-7)
  # test that conv_similarity is symmetric
  expect_equal(
    conv_similarity(midata, 1, 4, 1, cosine_sim),
    conv_similarity(midata, 4, 1, 1, cosine_sim))
  # for a vs b, we get the similarity a*a vs b
  expect_equal(conv_similarity(midata, 1, 2, 1, cosine_sim), 0.9889838)
})
