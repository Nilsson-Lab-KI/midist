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

test_that("filter_enrichment works correctly", {
  expect_equal(
    filter_enrichment(c(0.99, 0.01), tol = 0.01),
    c(0, 0))
  expect_equal(
    filter_enrichment(c(0.99, 0.01), tol = 0.01-1e-10),
    c(0.99, 0.01))
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

test_that("13-c correction sums to 1", {
  # 
  mid <- c(1,0,0)
  expect_equal(c13correct(mid), mid)
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


# TODO
#test_that("atom index is created correctly", {
#  expect_equal(
#    create_atom_index(c(1,1,3,3,3,5,6)),
#    
#  )
#})
 
 