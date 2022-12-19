
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
  Metabolite = c(rep("a",2), rep("b",3), rep("c",3), rep("d",4), rep("e",6)),
  Formula = c(rep("a",2), rep("b",3), rep("c",3), rep("d",4), rep("e",6)),
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
  expect_equal(dim(conv_reduce_all(midata, 1, cosine_dist, min, get_middle_met_matrix = F)), c(5,5))
  
  # make sure it's symmetric
  expect_true(isSymmetric(conv_reduce_all(midata, 1, cosine_dist, min, get_middle_met_matrix = F)))
  
  # for c (a zero vector) all cosine distances should be NAs due to division by zero
  expect_true(is.na(conv_reduce_all(midata, 1, cosine_dist, min, get_middle_met_matrix = F)[,3]))
  expect_true(is.na(conv_reduce_all(midata, 1, cosine_dist, min, get_middle_met_matrix = F)[3,]))
  
  expect_equal(as.numeric(format(conv_reduce_all(midata, 1, cosine_dist, min, get_middle_met_matrix = F)[2,1], digits = 4)),
               as.numeric(format(1.101622e-02, digits = 4)))
    
  # for a vs d, we take maximum of a*b vs d and a*c vs d
  expect_equal(
    conv_reduce_all(midata, 1, 4, 1, cosine_sim), 0.9949405, tolerance = 1e-7)
  # test that conv_similarity is symmetric
  expect_equal(
    conv_similarity(midata, 1, 4, 1, cosine_sim),
    conv_similarity(midata, 4, 1, 1, cosine_sim))
  # for a vs b, we get the similarity a*a vs b
  expect_equal(
    conv_similarity(midata, 1, 2, 1, cosine_sim), 0.9889838, tolerance = 1e-7)
  # for a vs c, we get a*a vs c which is NA
  expect_true(is.na(conv_similarity(midata, 1, 3, 1, cosine_sim)))
  # for a vs. e there is no 4-carbon peak to convolute with, return zero
  expect_equal(conv_similarity(midata, 1, 5, 1, cosine_sim), 0)
})


# test_that("conv_similarity returns a valid similarity score", {
#   # create MIData object
#   midata <- MIData(peak_areas_example, exp_names = "exp1")
#   # for b vs c (equal carbons) we get cosine_sim(b, c) which is NA
#   expect_true(is.na(conv_similarity(midata, 2, 3, 1, cosine_sim)))
#   # for a vs d, we take maximum of a*b vs d and a*c vs d
#   expect_equal(
#     conv_similarity(midata, 1, 4, 1, cosine_sim), 0.9949405, tolerance = 1e-7)
#   # test that conv_similarity is symmetric
#   expect_equal(
#     conv_similarity(midata, 1, 4, 1, cosine_sim),
#     conv_similarity(midata, 4, 1, 1, cosine_sim))
#   # for a vs b, we get the similarity a*a vs b
#   expect_equal(
#     conv_similarity(midata, 1, 2, 1, cosine_sim), 0.9889838, tolerance = 1e-7)
#   # for a vs c, we get a*a vs c which is NA
#   expect_true(is.na(conv_similarity(midata, 1, 3, 1, cosine_sim)))
#   # for a vs. e there is no 4-carbon peak to convolute with, return zero
#   expect_equal(conv_similarity(midata, 1, 5, 1, cosine_sim), 0)
# })
# 
# test_that("conv_distance returns a valid distance", {
#   # create MIData object
#   midata <- MIData(peak_areas_example, exp_names = "exp1")
#   # for b vs c (equal carbons) we get euclidean_dist(b, c)
#   expect_equal(
#     conv_distance(midata, 2, 3, 1, euclidean_dist), 0.8124038, tolerance = 1e-7)
#   # for a vs d, we take minimum of a*b vs d and a*c vs d
#   expect_equal(
#     conv_distance(midata, 1, 4, 1, euclidean_dist), 0.08158431, tolerance = 1e-7)
#   # test that conv_distance is symmetric
#   expect_equal(
#     conv_distance(midata, 1, 4, 1, euclidean_dist),
#     conv_distance(midata, 4, 1, 1, euclidean_dist))
# })
# 
# # example peak area data for testing conv_reduce_all, one experiment
# # here atom numbers are not increasing over peaks
# peak_areas_example2 <- data.frame(
#   peak_id <- c(rep("a",2), rep("b",3), rep("c",3), rep("d",4), rep("e",6)),
#   exp1 <- c(
#     0.8, 0.05, 0.1, 0.05,
#     0.8, 0.1, 0.1,
#     0.98, 0.02,
#     0.1, 0.0, 0.3, 0.0, 0.2, 0.05,
#     0.0, 0.0, 0.0)
# )
# 
# test_that("conv_reduce_all gives same values as conv_reduce", {
#   # create MIData object
#   midata <- MIData(peak_areas_example2, exp_names = "exp1")
#   # distance matrix from conv_reduce_all
#   dist_mat <- conv_reduce_all(midata, 1, euclidean_dist, min_dist)
#   # test row by row against conv_reduce
#   for(x in 1:5) {
#     dist_row <- sapply(1:5,
#                        function(y) conv_reduce(midata, x, y, 1, euclidean_dist, min_dist))
#     expect_equal(dist_mat[x,], dist_row)
#   }
# })
# 
