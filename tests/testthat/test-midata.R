#
# Tests for MIData objects
#

# example "peak area" data

peak_areas_1 <- data.frame(
  # two metabolites 'x' and 'y', 2 and 3 carbons respectively
  metabolite = c("x", "x", "x", "y", "y", "y", "y"),
  # we require a formula column
  annotation = "",
  # experiment 1 with 2 replicates
  exp1_1 = c(3, 2, 7, 9, 3, 0, 0),
  exp1_2 = c(4, 2, 6, 8, 3, 0, 0),
  # experiment 2 with 1 replicate
  exp2_1 = c(5, 0, 1, 2, 6, 5, 5)
)
# unique experiment names
exp_names_1 <- c('exp1', 'exp1', 'exp2')
# create MIData object
mi_data_1 <- MIData(peak_areas_1, exp_names_1)



test_that("MIData objects are created correctly", {
  # check object properties
  expect_equal(mi_data_1$peak_ids, c("x", "y"))
  expect_equal(length(mi_data_1$peak_index), 2)
  expect_equal(mi_data_1$peak_n_atoms, c(2, 3))
  expect_equal(mi_data_1$experiments, c("exp1", "exp2"))
  expect_equal(mi_data_1$exp_index, c(1, 3))
  expect_equal(mi_data_1$exp_n_rep, c(2, 1))
  expect_equal(mi_data_1$n_atoms_index[["2"]], c(1))
  expect_equal(mi_data_1$n_atoms_index[["3"]], c(2))

  # specifying columns
  mi_data_tmp <- MIData(peak_areas_1, exp_columns = 2 + c(1, 3))
  expect_equal(mi_data_tmp$experiments, c("exp1_1", "exp2_1"))
})


test_that("atom index is created correctly", {
  atom_index <- create_atom_index(c(1, 1, 3, 4, 4, 3, 6))
  expect_equal(
    atom_index,
    list("1" = c(1, 2), "3" = c(3, 6), "4" = c(4, 5), "6" = c(7))
  )
})


test_that("find_mi_index works correctly", {
  expect_equal(find_mi_index(c(3,2,1,5,2)), c(1,5,8,10,16))
})


test_that("get_exp_indices works correctly", {
  # create an MIData with multiple experiments
  mi_data <- MIData(
    data.frame(
      Metabolite = rep("a", 2),
      Formula = rep("a", 2),
      exp1 = c(0.80, 0.20),
      exp1 = c(0.10, 0.90),
      exp2 = c(0.15, 0.85),
      exp2 = c(0.50, 0.50),
      exp2 = c(0.40, 0.60),
      exp3 = c(0.10, 0.90),
      check.names = FALSE
    )
  )
  expect_equal(
    get_exp_indices(mi_data, 1),
    c(1, 2)
  )
  expect_equal(
    get_exp_indices(mi_data, 2),
    c(3, 4, 5)
  )
  expect_equal(
    get_exp_indices(mi_data, 3),
    c(6)
  )
})


test_that("get_peak_index works correctly", {
  expect_equal(get_peak_index(mi_data_1, "x"), 1)
  expect_equal(get_peak_index(mi_data_1, "y"), 2)
})


test_that("get_mids works correctly", {
  expect_equal(
    get_mids(mi_data_1, 1, 1),
    matrix(
      c(
        0.2500000, 0.3333333,
        0.1666667, 0.1666667,
        0.5833333, 0.5000000
      ),
      nrow = 3, byrow = TRUE
    ),
    tolerance = 1e-6
  )
  expect_equal(
    get_mids(mi_data_1, "x", 1),
    get_mids(mi_data_1, 1, 1)
  )

  expect_equal(
    get_mids(mi_data_1, 1, 2),
    matrix(
      c(
        0.8333333,
        0.0000000,
        0.1666667
      ),
      nrow = 3, byrow = TRUE
    ),
    tolerance = 1e-6
  )
  expect_equal(
    get_mids(mi_data_1, 1, "exp2"),
    get_mids(mi_data_1, 1, 2)
  )
  expect_equal(
    get_mids(mi_data_1, 2, 1),
    matrix(
      c(
        0.75, 0.7272727,
        0.25, 0.2727273,
        0.00, 0.0000000,
        0.00, 0.0000000
      ),
      nrow = 4, byrow = TRUE
    ),
    tolerance = 1e-6
  )
  expect_equal(
    get_mids(mi_data_1, "y", "exp1"),
    get_mids(mi_data_1, 2, 1)
  )
})


test_that("get_avg_mid works correctly", {
  # get a single MID vector
  expect_equal(
    get_avg_mid(mi_data_1, 1, 1),
    c(0.2916667, 0.1666667, 0.5416667),
    tolerance = 1e-6
  )
  expect_equal(
    get_avg_mid(mi_data_1, "x", "exp1"),
    get_avg_mid(mi_data_1, 1, 1)
  )
  expect_equal(
    get_avg_mid(mi_data_1, 1, 2),
    c(0.8333333, 0, 0.1666667),
    tolerance = 1e-6
  )
  expect_equal(
    get_avg_mid(mi_data_1, 1, 2),
    get_avg_mid(mi_data_1, "x", "exp2")
  )
  # specify list of experiments to get an MID matrix
  expect_equal(
    get_avg_mid(mi_data_1, 1, c(1, 2)),
    matrix(
      c(
        0.2916667, 0.8333333,
        0.1666667, 0.0000000,
        0.5416667, 0.1666667
      ),
      nrow = 3, byrow = TRUE
    ),
    tolerance = 1e-6
  )
  expect_equal(
    get_avg_mid(mi_data_1, "x", c("exp1", "exp2")),
    get_avg_mid(mi_data_1, 1, c(1, 2))
  )
  # leave out the exp parameter to get an MID matrix
  expect_equal(
    get_avg_mid(mi_data_1, 1),
    matrix(
      c(
        0.2916667, 0.8333333,
        0.1666667, 0.0000000,
        0.5416667, 0.1666667
      ),
      nrow = 3, byrow = TRUE
    ),
    tolerance = 1e-6
  )
  expect_equal(
    get_avg_mid(mi_data_1, "x"),
    get_avg_mid(mi_data_1, 1)
  )
})


test_that("midata_transform works correctly", {
  # apply identity function
  new_midata <- midata_transform(mi_data_1, identity)
  expect_equal(new_midata, mi_data_1)
  # apply 13C correction
  new_midata <- midata_transform(mi_data_1, c13correct)
  # check corrected MID
  expect_equal(
    get_mids(new_midata, 1, 1)[,1],
    c13correct(get_mids(mi_data_1, 1, 1)[,1])
  )
  # check averaged corrected MID
  expect_equal(
    get_avg_mid(new_midata, 1, 1),
    c13correct(get_avg_mid(mi_data_1, 1, 1))
  )
})


# an MIData object with 5 metabolites
peak_areas_2 <- data.frame(
  metabolite = c(rep("a",4), rep("b",3), rep("c",2), rep("d",6), rep("e",3)),
  annotation = "",
  exp1 = c(
    0.8, 0.05, 0.1, 0.05,
    0.8, 0.1, 0.1,
    0.98, 0.02,
    0.1, 0.0, 0.3, 0.0, 0.2, 0.05,
    0.0, 0.0, 0.0)
)
exp_names_2 = c("exp1")
mi_data_2 <- MIData(peak_areas_2, exp_names_2)


test_that("zero peaks are handled properly", {
  # zero peaks in peak areas should become NAs in the midata object
  expect_equal(all(is.na(get_mids(mi_data_2, 5, 1))), TRUE)
  # avg mids of NA mids should be assigned binomial
  expect_equal(
    get_avg_mid(mi_data_2, 5, 1),
    stats::dbinom(0:2, 2, natural_13C_fraction)
  )
})


test_that("all average mids sum to 1", {
  for(p in 1:length(mi_data_2$peak_ids)) {
    expect_equal(
      sum(get_avg_mid(mi_data_2, p, 1)),
      1
    )
  }
})


test_that("midata_subset works correctly", {
  # take a subset
  sub_index <- c(1,3,2)
  mi_data_sub <- midata_subset(mi_data_2, sub_index)
  expect_equal(mi_data_sub$peak_ids, mi_data_2$peak_ids[sub_index])
  expect_equal(mi_data_sub$peak_n_atoms, mi_data_2$peak_n_atoms[sub_index])
  expect_equal(
    get_avg_mid(mi_data_sub, 3, 1),
    get_avg_mid(mi_data_2, 2, 1))
  expect_equal(
    midata_subset(mi_data_2, c("a", "c", "b")),
    mi_data_sub
  )

  # permutation (all peaks, but in different order)
  sub_index <- c(4,2,1,5,3)
  mi_data_sub <- midata_subset(mi_data_2, sub_index)
  expect_equal(mi_data_sub$peak_ids, mi_data_2$peak_ids[sub_index])
  expect_equal(mi_data_sub$peak_n_atoms, mi_data_2$peak_n_atoms[sub_index])
  expect_equal(
    get_avg_mid(mi_data_sub, 4, 1),
    get_avg_mid(mi_data_2, 5, 1))
  expect_equal(
    midata_subset(mi_data_2, c("d", "b", "a", "e", "c")),
    mi_data_sub
  )

  # subset to entire range should give identical object
  expect_equal(midata_subset(mi_data_2, 1:5), mi_data_2)

})


# an MIData object with 3 metabolites across two experiments
peak_areas_3 <- data.frame(
  Metabolite = c(rep("a",4), rep("b",3), rep("c",3)),
  Formula = c(rep("a",4), rep("b",3), rep("c",3)),
  exp1 = c(
    0.8, 0.05, 0.1, 0.05,
    0.8, 0.1, 0.1,
    0.1, 0.4, 0.5),
  exp2 = c(
    0.85, 0.1, 0.0, 0.05,
    0.75, 0.15, 0.1,
    0.15, 0.35, 0.5)
)
exp_names_3 = c("exp1", "exp2")
mi_data_3 <- MIData(peak_areas_3, exp_names_3)


test_that("get_avg_mids works correctly", {
  # two peaks, single experiment gives an MI x peak matrix
  expect_equal(
    get_avg_mids(mi_data_3, c(2, 3), 1),
    matrix(
      c(
        0.8, 0.1, 0.1,
        0.1, 0.4, 0.5
      ),
      nrow = 3, byrow = FALSE
    ),
    tolerance = 1e-6
  )
  expect_equal(
    get_avg_mids(mi_data_3, c("b", "c"), "exp1"),
    get_avg_mids(mi_data_3, c(2, 3), 1)
  )
  expect_equal(
    get_avg_mids(mi_data_3, c(2, 3), 2),
    matrix(
      c(
        0.75, 0.15, 0.1,
        0.15, 0.35, 0.5
      ),
      nrow = 3, byrow = FALSE
    )
  )
  # two peaks, two experiments gives an MI x experiments x peaks array
  mid_array <- get_avg_mids(mi_data_3, c(2, 3), c(1, 2))
  expect_equal(
    dim(mid_array),
    c(3, 2, 2)
  )
  expect_equal(
    mid_array[, , 1],
    get_avg_mid(mi_data_3, 2, c(1, 2))
  )
  expect_equal(
    mid_array[, , 2],
    get_avg_mid(mi_data_3, 3, c(1, 2))
  )
  expect_equal(
    get_avg_mids(mi_data_3, c("b", "c"), c("exp1", "exp2")),
    mid_array
  )
})


#
# MID normalization -- move this to test-mid.R ?
#

test_that("normalized MIDs sum to 1", {
  # example "peak area" data, each column an MID
  peak_areas <- matrix(
    c(
      1, 0,
      3, 2,
      7, 0
    ),
    nrow = 3, ncol = 2
  )
  mids <- normalize_mids(peak_areas)
  expect_equal(colSums(mids), c(1, 1))
})

test_that("normalizing NA vector gives NA vector", {
  # example "peak area" data, each column an MID
  mids <- matrix(
    c(
      NA,
      NA,
      NA
    ),
    nrow = 3, ncol = 1
  )
  normalized_mids <- normalize_mids(mids)
  expect_equal(is.na(colSums(normalized_mids)), T)
})


test_that("false isotope removal works correctly", {
  mi_data_censored <- censor_false_mi(mi_data_1, threshold = 0.03, min_experiments = 2)
  # for peak 1 this removes M+2
  expect_equal(
    get_mids(mi_data_censored, 1, 1),
    matrix(
      c(
        0.6, 0.666667,
        0.4, 0.333333,
        0.0, 0.0
      ),
      nrow = 3, byrow = TRUE
    ),
    tolerance = 1e-6
  )
  expect_equal(
    get_avg_mid(mi_data_censored, 1),
    matrix(
      c(
        0.633333, 1.0,
        0.366667, 0.0,
        0.0, 0.0
      ),
      nrow = 3, byrow = TRUE
    ),
    tolerance = 1e-6
  )

  # for peak 2 this removes M+1
  expect_equal(
    get_avg_mid(mi_data_censored, 2),
    matrix(
      c(
        1.000, 0.166667,
        0.000, 0.000,
        0.000, 0.416667,
        0.000, 0.416667
      ),
      nrow = 4, byrow = TRUE
    ),
    tolerance = 1e-6
  )
})


test_that("misplace_peak_ids is correct", {
  # for this MIData object the only possible misplacement
  # is to exchange e and b
  expect_equal(
    misplace_peak_ids(mi_data_2),
    c("a", "e", "c", "d", "b")
  )
  # an MIData object with 20 peaks of atom sizes 1 and 2, 10 peaks of each
  peak_ids <- c(
    unlist(lapply(1:10, function(i) rep(paste0("x", i), 2))),
    unlist(lapply(1:10, function(i) rep(paste0("y", i), 3)))
  )
  unique_peak_ids <- unique(peak_ids)
  mi_data <- MIData(
    data.frame(
      Metabolite = peak_ids,
      Formula = peak_ids,
      exp1 = runif(2*10 + 3*10, 0, 1)
    )
  )
  misplaced_ids <- misplace_peak_ids(mi_data)
  # number of IDs must be the same
  expect_equal(length(misplaced_ids), 20)
  # misplacement occurs only within each atom size
  expect_equal(
    sort(misplaced_ids[1:10]),
    sort(unique_peak_ids[1:10])
  )
  expect_equal(
    sort(misplaced_ids[11:20]),
    sort(unique_peak_ids[11:20])
  )
  # each peak id differs from the original
  expect_false(any(misplaced_ids[1:10] == unique_peak_ids[1:10]))
  expect_false(any(misplaced_ids[11:20] == unique_peak_ids[11:20]))
})

