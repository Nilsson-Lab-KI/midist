# turn off warnings
options(warn = -1)

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
test_that("euclidean distance is correct", {
  # identical vectors have similarity 0
  expect_equal(euclidean_dist(c(2,3), c(2,3)), 0)
  # squared distance between unit vectors is 2
  expect_equal(euclidean_dist(c(1,0), c(0,1)), sqrt(2))
})
test_that("squared euclidean distance is correct", {
  # identical vectors have similarity 0
  expect_equal(euclidean_dist_sq(c(2,3), c(2,3)), 0)
  # squared distance between unit vectors is 2
  expect_equal(euclidean_dist_sq(c(1,0), c(0,1)), 2)
})


# example peak area data for testing conv_reduce_all, one experiment
# example-1: carbon numbers are sorted
peak_areas_example <- data.frame(
  Metabolite = c(rep("a",2), rep("b",3), rep("c",3), rep("d",4), rep("e",6)),
  Formula = c(rep("a",2), rep("b",3), rep("c",3), rep("d",4), rep("e",6)),
  exp1 = c(
    0.98, 0.02,
    0.8, 0.1, 0.1,
    0.0, 0.0, 0.0, # zero norm case
    0.8, 0.05, 0.1, 0.05,
    0.1, 0.0, 0.3, 0.0, 0.2, 0.05))


# example-2: here atom numbers are not increasing over peaks
peak_areas_example2 <- data.frame(
  Metabolite = c(rep("a",4), rep("b",3), rep("c",2), rep("d",6), rep("e",3)),
  Formula = c(rep("a",4), rep("b",3), rep("c",2), rep("d",6), rep("e",3)),
  exp1 = c(
    0.8, 0.05, 0.1, 0.05,
    0.8, 0.1, 0.1,
    0.98, 0.02,
    0.1, 0.0, 0.3, 0.0, 0.2, 0.05,
    0.0, 0.0, 0.0)
)

test_that("conv_reduce_all returns a valid distance for example-1", {
  # create MIData object
  midata <- MIData(peak_areas_example, exp_names = "exp1")
  
  # make sure it outputs a matrix of correct size
  expect_equal(dim(conv_reduce_all(midata, 1, cosine_sim, max, get_middle_met_matrix = F)), c(5,5))
  
  # make sure it's symmetric
  expect_true(isSymmetric(conv_reduce_all(midata, 1, cosine_sim, max, get_middle_met_matrix = F)))
  
  # for c (a zero vector) all cosine distances should be NAs due to division by zero
  expect_true(unique(is.na(conv_reduce_all(midata, 1, cosine_sim, max, get_middle_met_matrix = F, which.max)[,3])))
  expect_true(unique(is.na(conv_reduce_all(midata, 1, cosine_sim, max, get_middle_met_matrix = F, which.max)[3,])))
  
  expect_equal(as.numeric(format(conv_reduce_all(midata, 1, cosine_dist, min, get_middle_met_matrix = F)[2,1], digits = 4)),
               as.numeric(format(1.101622e-02, digits = 4)))
    
  # for a vs d, we take maximum of a*b vs d and a*c vs d
  expect_equal(
    conv_reduce_all(midata, 1, cosine_sim, max, F, which.max)[4,2], 0.9949405, tolerance = 1e-7)
  # test that conv_reduce_all returns a symmetric matrix
  expect_equal(
    conv_reduce_all(midata, 1, cosine_sim, max, F, which.max),
    t(conv_reduce_all(midata, 1, cosine_sim, max, F, which.max)))
  # for a vs b, we get the similarity a*a vs b
  expect_equal(
    conv_reduce_all(midata, 1, cosine_sim, max, F, which.max)[2,1], 0.9889838, tolerance = 1e-7)
  # for a vs c, we get a*a vs c which is NA
  expect_true(is.na(conv_reduce_all(midata, 1, cosine_sim, max, F, which.max)[3,1]))
  # for a vs. e there is no 4-carbon peak to convolute with, return inf
  expect_true(is.infinite(conv_reduce_all(midata, 1, cosine_sim, max, F, which.max)[5,1]))
})
test_that("conv_reduce_all returns a valid middle metabolite matrix for example-1", {
  # create MIData object
  midata <- MIData(peak_areas_example, exp_names = "exp1")
  # compute similarity and middle metabolite matrices
  output <- conv_reduce_all(midata, 1, cosine_sim, max, get_middle_met_matrix = T, which.max)
  
  # make sure the output is a list
  expect_true(is.list(output))
  # of length 2 (one for similarity and one for the middle metabolite)
  expect_equal(length(output), 2)
  
  # C pool: available C groups
  c_pool <- midata$peak_n_atoms
  # C number difference for all pairs as a matrix of the same size as the output matrices
  c_diff <- matrix(abs(expand.grid(midata$peak_n_atoms, midata$peak_n_atoms)[,1] - expand.grid(midata$peak_n_atoms, midata$peak_n_atoms)[,2]),
                   5, 5)
  # expect diagonals to be zero
  expect_equal(unique(diag(c_diff)), 0)
  # get the index of non non zero values and those that we have in our carbon pool
  index <- which(c_diff != 0 & c_diff %in% c_pool)
  # make sure we only made predictions for these metabolites
  expect_equal(length(setdiff(which(is.na(output[[2]]) == F), index)), 0)
})
test_that("conv_reduce_all returns a valid distance for example-2", {
  # create MIData object
  midata <- MIData(peak_areas_example2, exp_names = "exp1")
  
  # make sure it outputs a matrix of correct size
  expect_equal(dim(conv_reduce_all(midata, 1, cosine_dist, min, get_middle_met_matrix = F)), c(5,5))
  
  # make sure it's symmetric
  expect_true(isSymmetric(conv_reduce_all(midata, 1, cosine_dist, min, get_middle_met_matrix = F)))
  
  # for e (a zero vector) all cosine distances should be NAs due to division by zero
  expect_true(unique(is.na(conv_reduce_all(midata, 1, cosine_dist, min, get_middle_met_matrix = F)[,5])))
  expect_true(unique(is.na(conv_reduce_all(midata, 1, cosine_dist, min, get_middle_met_matrix = F)[5,])))
  
  # expect a valid distance
  expect_equal(as.numeric(sprintf(conv_reduce_all(midata, 1, cosine_dist, min, get_middle_met_matrix = F)[2,1], fmt = '%#.5f')),
               0.00506)
  
  # for c vs a, we take maximum of c*b vs a and c*e vs a
  expect_equal(
    conv_reduce_all(midata, 1, cosine_sim, max, F, which.max)[4,2], 0.4587568, tolerance = 1e-7)
  # test that conv_reduce_all returns a symmetric matrix
  expect_equal(
    conv_reduce_all(midata, 1, cosine_sim, max, F, which.max),
    t(conv_reduce_all(midata, 1, cosine_sim, max, F, which.max)))
  # for a vs b, we get the similarity a*a vs b
  expect_equal(
    conv_reduce_all(midata, 1, cosine_sim, max, F, which.max)[2,1], 0.9949405, tolerance = 1e-7)
  # for a vs c, we get a*a vs c which is NA
  expect_true(is.na(conv_reduce_all(midata, 1, cosine_sim, max, F, which.max)[3,1]))
  # for c vs. d there is no 4-carbon peak to convolute with, return inf
  expect_true(is.infinite(conv_reduce_all(midata, 1, cosine_sim, max, F, which.max)[3,4]))
})
test_that("conv_reduce_all returns a valid middle metabolite matrix for example-2", {
  # create MIData object
  midata <- MIData(peak_areas_example2, exp_names = "exp1")
  # compute similarity and middle metabolite matrices
  output <- conv_reduce_all(midata, 1, cosine_sim, max, get_middle_met_matrix = T, which.max)
  
  # make sure the output is a list
  expect_true(is.list(output))
  # of length 2 (one for similarity and one for the middle metabolite)
  expect_equal(length(output), 2)
  
  # C pool: available C groups
  c_pool <- midata$peak_n_atoms
  # C number difference for all pairs as a matrix of the same size as the output matrices
  c_diff <- matrix(abs(expand.grid(midata$peak_n_atoms, midata$peak_n_atoms)[,1] - expand.grid(midata$peak_n_atoms, midata$peak_n_atoms)[,2]),
                   5, 5)
  # expect diagonals to be zero
  expect_equal(unique(diag(c_diff)), 0)
  # get the index of non non zero values and those that we have in our carbon pool
  index <- which(c_diff != 0 & c_diff %in% c_pool)
  # make sure we only made predictions for these metabolites
  expect_equal(length(setdiff(which(is.na(output[[2]]) == F), index)), 0)
})
test_that("conv_reduce_all gives same values as conv_reduce", {
  # create MIData object
  midata <- MIData(peak_areas_example2, exp_names = "exp1")
  # distance matrix from conv_reduce_all
  dist_mat <- conv_reduce_all(midata, 1, euclidean_dist, min, F, which.min)
  # test row by row against conv_reduce
  for(x in 1:5) {
    dist_row <- sapply(1:5,
                       function(y) conv_reduce(midata, x, y, 1, euclidean_dist, min))
    expect_equal(dist_mat[x,], dist_row)
  }
})

