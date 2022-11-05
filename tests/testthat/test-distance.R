
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
  #
  sim_data <- list(peak_ids = c("a", "b", "c", "d"),
                   peak_index = c(1, 3, 6, 9),
                   peak_n_atoms = c(1, 2, 2, 3),
                   n_atoms_index = list('1' = 1, '2' = c(2, 3), '3' = 4),
                   experiments = c(1,2), exp_n_rep = c(1,1),
                   peak_areas = c(0.98, 0.02, 0.8, 0.1, 0.1, 0.6, 0.2, 0,2, 0.8, 0.05, 0.1, 0.05),
                   mids = cbind(c(0.98, 0.02, 0.8, 0.1, 0.1, 0.6, 0.2, 0,2, 0.8, 0.05, 0.1, 0.05),
                                c(0.98, 0.02, 0.8, 0.1, 0.1, 0.6, 0.2, 0,2, 0.8, 0.05, 0.1, 0.05)),
                   avg_mids = cbind(c(0.98, 0.02, 0.8, 0.1, 0.1, 0.6, 0.2, 0,2, 0.8, 0.05, 0.1, 0.05),
                                    c(0.98, 0.02, 0.8, 0.1, 0.1, 0.6, 0.2, 0,2, 0.8, 0.05, 0.1, 0.05)))

  expect_equal(round(conv_similarity(sim_data, 1, 4, 1, cosine_sim), 7), 0.9979649)
})
