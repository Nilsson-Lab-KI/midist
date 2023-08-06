
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

# test conv_reduce on example 1 with given f and g functions
test_conv_reduce_1 <- function(f, g)
{
  # compute full matrix row by row
  # NOTE: assignment to matrix element discards attributes
  conv_mat <- matrix(NA, nrow = 5, ncol = 5)
  for(x in 1:5) {
    for(y in 1:5) {
      conv_mat[x,y] <- conv_reduce(midata_1, x, y, 1, f, g)
    }
  }

  # make sure matrix is symmetric
  expect_true(isSymmetric(conv_mat))

  # distance for a vs b = f(a*a, b), selects a
  expect_equal(conv_reduce(midata_1, 1, 2, 1, f, g),
               g(with_attr(
                 c(f(convolute(a_mid, a_mid), b_mid)), "index", c(1))),
               tolerance = 1e-7)
  # distance for a vs c = f(a*a, c), where c is a zero vector, selects a
  expect_equal(conv_reduce(midata_1, 1, 3, 1, f, g),
               g(with_attr(
                 c(f(convolute(a_mid, a_mid), c_mid)), "index", c(1))),
               tolerance = 1e-7)
  # distance for a vs d = g(f(a*b, d), f(a*c, d)), selects b
  expect_equal(conv_reduce(midata_1, 1, 4, 1, f, g),
               g(with_attr(
                 c(
                   f(d_mid, convolute(a_mid, b_mid)),
                   f(d_mid, convolute(a_mid, c_mid))), "index", c(2,3))),
               tolerance = 1e-7)
  # for a vs. e there is no 4-carbon peak to convolute with,
  # so expect g(c()), index = NA
  expect_equal(conv_reduce(midata_1, 1, 5, 1, f, g),
               g(c()),
               tolerance = 1e-7)
  # distance for b vs d = f(a*b, d), select a
  expect_equal(conv_reduce(midata_1, 2, 4, 1, f, g),
               g(with_attr(
                 c(f(convolute(a_mid, b_mid), d_mid)), "index", c(1))),
               tolerance = 1e-7)
  return(conv_mat)
}


test_that("conv_reduce is correct on minimum euclidean distance", {
  test_conv_reduce_1(euclidean_dist, min_nonempty)
})


test_that("conv_reduce is correct on max cosine similarity", {
  # maximum cosine similarity
  # in this case, the zero vector c causes NA values which must be handled
  sim_mat <- test_conv_reduce_1(cosine_sim, max_nonempty)
  # for c (a zero vector) all cosine distances should be NAs due to division by zero
  expect_true(all(is.na(sim_mat[,3])))
  expect_true(all(is.na(sim_mat[3,])))
})

test_that("conv_reduce_all returns a matrix with index attribute", {
  # maximum cosine similarity
  conv_all_mat <- conv_reduce_all(midata_1, 1, cosine_sim, max_nonempty)
  # value matrix
  expect_true(is.matrix(conv_all_mat))
  expect_equal(dim(conv_all_mat), c(5, 5))
  # index matrix
  expect_true(is.matrix(attr(conv_all_mat, "index")))
  expect_equal(dim(attr(conv_all_mat, "index")), c(5, 5))
})


# compare conv_reduce_all elementwise to conv_reduce on the given mi_data
# with functions f and g
test_conv_reduce_all <- function(mi_data, f, g)
{
  # matrix from conv_reduce_all
  conv_all_mat <- conv_reduce_all(mi_data, 1, f, g)
  # make sure matrix is symmetric
  expect_true(isSymmetric(conv_all_mat))

  # check against conv_reduce element by element
  n_peaks <- length(mi_data$peak_ids)
  for(x in 1:n_peaks) {
    for(y in 1:n_peaks) {
      conv_all_elem <- with_attr(conv_all_mat[x,y], "index",
                                 attr(conv_all_mat, "index")[x,y])
      expect_equal(conv_all_elem,
                   conv_reduce(mi_data, !!x, !!y, 1, f, g))
    }
  }
}


test_that("conv_reduce_all returns a valid similarity for example-1", {
  # maximum cosine similarity
  test_conv_reduce_all(midata_1, cosine_sim, max_nonempty)
})


test_that("conv_reduce_all returns a valid distance for example-1", {
  test_conv_reduce_all(midata_1, euclidean_dist, min_nonempty)
})


# example-2:  as example 1 but peaks are permuted
midata_2 <- midata_subset(midata_1, c(4,2,1,5,3))

test_that("conv_reduce_all returns a valid distance matrix for example-2", {
  test_conv_reduce_all(midata_2, cosine_dist, min_nonempty)
})


# test case from HMEC data
adp_mid <- c(0.3739468302, 0.0868470331, 0.1394170159, 0.1083851249,
             0.2638703387, 0.0266737617, 0.0008598955, 0.0000000000,
             0.0000000000, 0.0000000000, 0.0000000000)

gthrd_mid <- c(2.910509e-01, 4.279995e-02, 6.118703e-01, 5.421727e-02,
               4.515372e-06, 5.363846e-05, 0.000000e+00, 3.399253e-06,
               0.000000e+00, 0.000000e+00, 0.000000e+00)

ibutcrn_mid <- c(8.863661e-01, 1.040638e-01, 9.447574e-03, 0.000000e+00,
                 6.969751e-05, 7.760243e-06, 0.000000e+00, 0.000000e+00,
                 0.000000e+00, 0.000000e+00, 4.509453e-05, 0.000000e+00)

nad1_mid <- c(3.185700e-01, 1.165741e-01, 1.398212e-01, 1.157578e-01,
              2.549923e-01, 5.174851e-02, 2.379705e-03, 8.610149e-05,
              0.000000e+00, 7.030201e-05, 0.000000e+00, 0.000000e+00,
              0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
              0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
              0.000000e+00, 0.000000e+00)

nad2_mid <- c(3.116997e-01, 1.203842e-01, 1.414667e-01, 1.142038e-01,
              2.490341e-01, 5.413205e-02, 8.889942e-03, 9.395356e-05,
              0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
              0.000000e+00, 0.000000e+00, 0.000000e+00, 4.062245e-05,
              5.491461e-05, 0.000000e+00, 0.000000e+00, 0.000000e+00,
              0.000000e+00, 0.000000e+00)

# example-1: carbon numbers are sorted
peak_areas_example_hmec <- data.frame(
  Metabolite = c(rep("adp", 10+1), rep("gthrd" ,10+1), rep("ibutcrn", 11+1), rep("nad1", 21+1), rep("nad2", 21+1)),
  Formula = c(rep("adp", 10+1), rep("gthrd" ,10+1), rep("ibutcrn", 11+1), rep("nad1", 21+1), rep("nad2", 21+1)),
  exp1 = c(adp_mid, gthrd_mid, ibutcrn_mid, nad1_mid, nad2_mid)
)

# create MIData object
midata_hmec <- MIData(peak_areas_example_hmec, exp_names = "exp1")

test_that("conv_reduce_all returns a valid matrix for HMEC example", {
  test_conv_reduce_all(midata_hmec, euclidean_dist_sq, min_nonempty)
})


# combining distance matrices across experiments
test_that("combine works correcly", {
  # example distance matrices
  dm_1 <- matrix(
    c(
      0.0,  0.2,  0.9,  0.4,
      0.2,  0.0,  NA,   0.1,
      0.9,  NA,   0.0,  0.1,
      0.4,  0.1,  0.1,  0.0),
    nrow = 4
  )
  dm_2 <- matrix(
    c(
      0.0,  0.8,  0.6,  0.5,
      0.8,  0.0,  NA,   0.2,
      0.6,  NA,   0.0,  0.1,
      0.5,  0.2,  0.1,  0.0),
    nrow = 4
  )
  # corresponding middle metabolite matrices
  mm_1 <- matrix(
    c(
      NA,  NA,  1,   1,
      NA,  NA,  NA,  2,
      1,   NA,  NA,  4,
      1,   2,   4,   NA),
    nrow = 4
  )
  mm_2 <- matrix(
    c(
      NA,  NA,  2,   1,
      NA,  NA,  NA,  3,
      2,   NA,  NA,  4,
      1,   3,   4,   NA),
    nrow = 4
  )

  expect_equal(
    remn::combine(list(dm_1, dm_2), list(mm_1, mm_2), max),
    list(
      matrix(
        c(
          0.0,  0.8,  0.9,  0.5,
          0.8,  0.0,  NA,   0.2,
          0.9,  NA,   0.0,  0.1,
          0.5,  0.2,  0.1,  0.0),
        nrow = 4
      ),
      matrix(
        c(
          NA,  NA,  1,   1,
          NA,  NA,  NA,  3,
          1,   NA,  NA,  4,
          1,   3,   4,   NA),
        nrow = 4
      ),
      # the experiment index is not unique, we might allow
      # any index for positions where there are ties ...
      matrix(
        c(
          1,   2 ,  1,   2,
          2,   1,   NA,  2,
          1,   NA,  1,   1,
          2,   2,   1,   1),
        nrow = 4
      )
    )
  )
})

# An MIData object with two experiments
peak_areas_3 <- data.frame(
  Metabolite = c(rep("a",2), rep("b",3), rep("c",3), rep("d",4), rep("e",6)),
  Formula = c(rep("a",2), rep("b",3), rep("c",3), rep("d",4), rep("e",6)),
  exp1 = c(
    0.98, 0.02,
    0.8, 0.1, 0.1,
    0.0, 0.8, 0.2,
    0.8, 0.05, 0.1, 0.05,
    0.1, 0.0, 0.3, 0.0, 0.2, 0.0
  ),
  exp2 = c(
    0.91, 0.09,
    0.2, 0.1, 0.7,
    0.4, 0.6, 0.0,
    0.1, 0.05, 0.1, 0.75,
    0.2, 0.0, 0.2, 0.0, 0.0, 0.2
  )
)
mi_data_3 <- MIData(peak_areas_3)


test_that("enrichment_dist_matrix is correct", {

  ed_matrix <- enrichment_dist_matrix(mi_data_3, c("exp1", "exp2"), "euclidean")
  expect_true(is.matrix(ed_matrix))
  expect_equal(
    ed_matrix,
    matrix(
      c(
        0.0000000, 0.67268120, 0.6168468, 0.75192346, 0.5842849,
        0.6726812, 0.00000000, 0.6363961, 0.08498366, 0.4249183,
        0.6168468, 0.63639610, 0.0000000, 0.70867639, 0.2134375,
        0.7519235, 0.08498366, 0.7086764, 0.00000000, 0.4955356,
        0.5842849, 0.42491829, 0.2134375, 0.49553562, 0.0000000
      ),
      nrow = 5,
      dimnames = list(letters[1:5], letters[1:5])
    ),
    tolerance = 1e-6
  )
})


test_that("remn_v1 is correct", {
  # list of distance matrices
  dm_list <- remn_v1(
      mi_data_3, f = euclidean_dist, g = min_nonempty, rdata_fname = "")
  expect_true(is.list(dm_list))
  expect_equal(length(dm_list), 2)
  # check matrices and attributes
  expect_equal(
    dm_list[[1]],
    with_attr(
      matrix(
        c(
          0.00000000, 0.19835564, 1.2413803, 0.08158431, NA,
          0.19835564, 0.00000000, 1.0677078, 0.08158431, 0.6745781,
          1.24138026, 1.06770783, 0.0000000, 1.09243581, 0.7812099,
          0.08158431, 0.08158431, 1.0924358, 0.00000000, 0.6745781,
          NA,         0.67457806, 0.7812099, 0.67457806, 0.0000000
        ),
        nrow = 5
      ),
      "index",
      matrix(
        c(
          NA,  1,   1,   2,  NA,
          1,   NA, NA,   1,   4,
          1,   NA, NA,   1,   4,
          2,   1,   1,  NA,   2,
          NA,  4,   4,   2,  NA
        ),
        nrow = 5
      )
    ),
    tolerance = 1e-6
  )
  expect_equal(
    dm_list[[2]],
    with_attr(
      matrix(
        c(
          0.0000000, 0.9366460, 0.6112329, 0.8833403, NA,
          0.9366460, 0.0000000, 0.8831761, 0.8833403, 0.5011653,
          0.6112329, 0.8831761, 0.0000000, 0.9577766, 0.7778175,
          0.8833403, 0.8833403, 0.9577766, 0.0000000, 0.5011653,
          NA,        0.5011653, 0.7778175, 0.5011653, 0.0000000
        ),
        nrow = 5
      ),
      "index",
      matrix(
        c(
          NA,  1,   1,   2,  NA,
          1,   NA, NA,   1,   4,
          1,   NA, NA,   1,   4,
          2,   1,   1,  NA,   2,
          NA,  4,   4,   2,  NA
        ),
        nrow = 5
      )
    ),
    tolerance = 1e-6
  )

})


test_that("remn_v2 is correct", {
  # a list [distance matrix, middle metabolite matrix]
  dm_mm_list <- remn_v2(
    mi_data_3, f = euclidean_dist, g = min_nonempty, rdata_fname = "", return = T)
  # distance matrix
  dm = dm_mm_list$distance_matrix
  expect_equal(
    dm,
    matrix(
      c(
        NA,        1.1350016, 1.852613, 0.9649246, NA,
        1.1350016, NA,        1.950884, 0.9649246, 1.175743,
        1.8526132, 1.9508839, NA,       2.0502124, 1.559027,
        0.9649246, 0.9649246, 2.050212, NA,        1.175743,
        NA,        1.1757434, 1.559027, 1.1757434, NA
      ),
      nrow = 5, dimnames = list(mi_data_3$peak_ids,mi_data_3$peak_ids)
    ),
    tolerance = 1e-6
  )
  # distance matrix
  mm = dm_mm_list$middle_metabolite_matrix
  expect_equal(
    mm,
    matrix(
      c(
        NA, 1,  1,  2,  NA,
        1,  NA, NA, 1,  4,
        1,  NA, NA, 1,  4,
        2,  1,  1,  NA, 2,
        NA, 4,  4,  2,  NA
      ),
      nrow = 5, dimnames = list(mi_data_3$peak_ids,mi_data_3$peak_ids)
    )
  )
})


