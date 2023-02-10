

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

# individual MIDs
a_mid <- c(0.98, 0.02)
b_mid <- c(0.8, 0.1, 0.1)
c_mid <- c(0.0, 0.0, 0.0)           # zero vector
d_mid <- c(0.8, 0.05, 0.1, 0.05)
e_mid <- c(0.1, 0.0, 0.3, 0.0, 0.2, 0.05)

# example-1: carbon numbers are sorted
peak_areas_example_1 <- data.frame(
  Metabolite = c(rep("a",2), rep("b",3), rep("c",3), rep("d",4), rep("e",6)),
  Formula = c(rep("a",2), rep("b",3), rep("c",3), rep("d",4), rep("e",6)),
  exp1 = c(a_mid, b_mid, c_mid, d_mid, e_mid)
)

# create MIData object
midata_1 <- MIData(peak_areas_example_1, exp_names = "exp1")


test_that("conv_reduce is correct on euclidean distance", {
  # minimum euclidean distance
  f <- euclidean_dist
  g <- min
  # distance for b vs a = f(b, a*a)
  expect_equal(
    conv_reduce(midata_1, 2, 1, 1, f, g),
    f(b_mid, convolute(a_mid, a_mid)),
    tolerance = 1e-7)
  # distance for a vs d = max(f(d, a*b), f(d, a*c))
  expect_equal(
    conv_reduce(midata_1, 1, 4, 1, f, g),
    g(
      f(d_mid, convolute(a_mid, b_mid)),
      f(d_mid, convolute(a_mid, c_mid))),
    tolerance = 1e-7)
})


test_that("conv_reduce is correct on cosine similarity", {
  # maximum cosine similarity
  # in this case, the zero vector c causes NA values which must be handled
  f <- cosine_sim
  g <- max
  # distance for b vs a = f(b, a*a)
  expect_equal(
    conv_reduce(midata_1, 2, 1, 1, f, g),
    f(b_mid, convolute(a_mid, a_mid)),
    tolerance = 1e-7)
  # for a vs d = max(f(d, a*b), f(d, a*c)) = f(d, a*b)
  # since f(d, a*c) is NA
  expect_equal(
    conv_reduce(midata_1, 1, 4, 1, f, g),
    f(d_mid, convolute(a_mid, b_mid)),
    tolerance = 1e-7)
})


test_that("conv_reduce_all returns a valid similarity for example-1", {
  # maximum cosine similarity
  # in this case, the zero vector c causes NA values which must be handled
  f <- cosine_sim
  g <- max

  sm_max <- conv_reduce_all(midata_1, 1, f, g, get_middle_met_matrix = F)
  # make sure matrix has correct size
  expect_equal(dim(sm_max), c(5,5))
  # make sure matrix is symmetric
  expect_true(isSymmetric(sm_max))
  # for c (a zero vector) all cosine distances should be NAs due to division by zero
  expect_true(all(is.na(sm_max[,3])))
  expect_true(all(is.na(sm_max[3,])))

  # distance for a vs d = min(f(d, a*b), f(d, a*c)) = f(d, a*b)
  # since f(d, a*c) is NA
  # TODO: this fails because conv_select all does not handle NA values correctly
  expect_equal(
    sm_max[1,4],
    f(d_mid, convolute(a_mid, b_mid)),
    tolerance = 1e-7)

  # for b vs d = f(d, a*b)
  expect_equal(
    sm_max[4,2],
    f(d_mid, convolute(a_mid, b_mid)),
    tolerance = 1e-7)
  # for a vs b = f(a*a, b)
  expect_equal(
    sm_max[2,1],
    f(convolute(a_mid, a_mid), b_mid),
    tolerance = 1e-7)
  # for a vs c, we get a*a vs c which is NA
  expect_true(is.na(sm_max[3,1]))
  # for a vs. e there is no 4-carbon peak to convolute with, return inf
  expect_true(is.infinite(sm_max[5,1]))
})


test_that("conv_reduce_all returns a valid distance for example-1", {
  # compute distance matrix with g = min
  dm_min <- conv_reduce_all(midata_1, 1, cosine_dist, min, get_middle_met_matrix = F)
  # make sure matrix is symmetric
  expect_true(isSymmetric(dm_min))
  # distance for b vs a = 0.01101622
  expect_equal(
    dm_min[2,1],
    cosine_dist(b_mid, convolute(a_mid, a_mid)),
    tolerance = 1e-7)
})


test_that("conv_reduce_all returns a valid middle metabolite matrix for example-1", {
  # compute similarity and middle metabolite matrices
  f <- cosine_sim
  g <- max
  output <- conv_reduce_all(midata_1, 1, f, g, get_middle_met_matrix = T, which.max)

  # make sure the output is a list
  expect_true(is.list(output))
  # of length 2 (one for similarity and one for the middle metabolite)
  expect_equal(length(output), 2)
  #print(output[[2]])
  n_atoms <- midata_1$peak_n_atoms
  # atom number difference for all pairs of metabolites
  # this is a matrix of the same size as the output matrices
  c_diff <- abs(outer(n_atoms, n_atoms, "-"))
  # get the index of nonzero values and those that we have in our carbon pool
  index <- which(c_diff != 0 & c_diff %in% n_atoms)
  # make sure we only made predictions for these metabolites
  expect_equal(
    length(setdiff(which(is.na(output[[2]]) == F), index)),
    0)
  # TODO: test that selected metabolite is correct
})


# example-2: here atom numbers are not increasing over peaks

# individual MIDs
a_mid <- c(0.8, 0.05, 0.1, 0.05)
b_mid <- c(0.8, 0.1, 0.1)
c_mid <- c(0.98, 0.02)
d_mid <- c(0.1, 0.0, 0.3, 0.0, 0.2, 0.05)
e_mid <- c(0.0, 0.0, 0.0)           # zero vector

peak_areas_example2 <- data.frame(
  # metabolite and #atoms: a 3, b 2, c 1, d 5, e 2
  Metabolite = c(rep("a",4), rep("b",3), rep("c",2), rep("d",6), rep("e",3)),
  # RN: not sure what Formula is used for here?
  Formula = c(rep("a",4), rep("b",3), rep("c",2), rep("d",6), rep("e",3)),
  exp1 = c(a_mid, b_mid, c_mid, d_mid, e_mid)
)

midata_2 <- MIData(peak_areas_example2, exp_names = "exp1")


test_that("conv_reduce_all returns a valid distance matrix for example-2", {
  # compute distance matrix for experiment 1 with g = min
  dm_min <- conv_reduce_all(midata_2, 1, f = cosine_dist, g = min, get_middle_met_matrix = F)

  # make sure distance matrix of correct size
  expect_equal(dim(dm_min), c(5,5))
  # make sure the distance matrix is symmetric
  expect_true(isSymmetric(dm_min))
  # for MID e (a zero vector) all cosine distances should be NAs due to division by zero
  expect_true(all(is.na(dm_min[,5])))
  expect_true(all(is.na(dm_min[5,])))

  # distance for b vs a = dist(a, b * c)
  expect_equal(dm_min[2,1],
               cosine_dist(a_mid, convolute(b_mid, c_mid)),
               tolerance = 1e-7)

  # compute similarity matrix with g = max
  sm_max <- conv_reduce_all(midata_2, 1, cosine_sim, max, F, which.max)
  # similarity for d vs b = sim(d, a * b)
  expect_equal(sm_max[4,2],
               cosine_sim(d_mid, convolute(a_mid, b_mid)),
               tolerance = 1e-7)

  # test that conv_reduce_all returns a symmetric matrix
  expect_true(isSymmetric(sm_max))
  # for a vs b, we get the similarity a*a vs b
  expect_equal(sm_max[2,1], 0.9949405, tolerance = 1e-7)
  # for a vs c, we get a*a vs c which is NA
  expect_true(is.na(sm_max[3,1]))
  # for c vs. d there is no 4-carbon peak to convolute with, return inf
  expect_true(is.infinite(sm_max[3,4]))
})


test_that("conv_reduce_all returns a valid middle metabolite matrix for example-2", {
  # compute similarity and middle metabolite matrices
  output <- conv_reduce_all(midata_2, 1, cosine_sim, max, get_middle_met_matrix = T, which.max)

  # make sure the output is a list
  expect_true(is.list(output))
  # of length 2 (one for similarity and one for the middle metabolite)
  expect_equal(length(output), 2)

  # C pool: available C groups
  c_pool <- midata_2$peak_n_atoms
  # C number difference for all pairs as a matrix of the same size as the output matrices
  c_diff <- matrix(
    abs(
      expand.grid(midata_2$peak_n_atoms, midata_2$peak_n_atoms)[,1]
      - expand.grid(midata_2$peak_n_atoms, midata_2$peak_n_atoms)[,2]),
    5, 5)
  # expect diagonals to be zero
  expect_equal(unique(diag(c_diff)), 0)
  # get the index of non non zero values and those that we have in our carbon pool
  index <- which(c_diff != 0 & c_diff %in% c_pool)
  # make sure we only made predictions for these metabolites
  expect_equal(length(setdiff(which(is.na(output[[2]]) == F), index)), 0)
})


test_that("conv_reduce_all gives same values as conv_reduce", {
  # distance matrix from conv_reduce_all
  dist_mat <- conv_reduce_all(midata_2, 1, euclidean_dist, min, F, which.min)
  # test row by row against conv_reduce
  for(x in 1:5) {
    dist_row <- sapply(1:5,
                       function(y) conv_reduce(midata_2, x, y, 1, euclidean_dist, min))
    expect_equal(dist_mat[x,], dist_row)
  }
})


