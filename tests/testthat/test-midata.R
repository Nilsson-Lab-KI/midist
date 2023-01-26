#
# Tests for MIData objects
#

# example "peak area" data
peak_areas <- data.frame(
  # two metabolites 'x' and 'y', 2 and 3 carbons respectively
  Metabolite = c("x", "x", "x", "y", "y", "y", "y"),
  # we require a formula column
  Formula = c("x", "x", "x", "y", "y", "y", "y"),
  # experiment 1 with 2 replicates
  exp1_1 = c(3, 2, 7, 9, 3, 0, 0),
  exp1_2 = c(4, 2, 6, 8, 3, 0, 0),
  exp2_1 = c(5, 0, 1, 2, 6, 5, 5)
)
# unique experiment names
exp_names <- c("exp1", "exp1", "exp2")


test_that("MIData objects are created correctly", {
  # create MIData object
  midata <- MIData(peak_areas, exp_names)
  # check object properties
  expect_equal(midata$peak_ids, c("x", "y"))
  expect_equal(length(midata$peak_index), 2)
  expect_equal(midata$peak_n_atoms, c(2, 3))
  expect_equal(midata$experiments, c("exp1", "exp2"))
  expect_equal(midata$exp_index, c(1, 3))
  expect_equal(midata$exp_n_rep, c(2, 1))
  expect_equal(midata$n_atoms_index[["2"]], c(1))
  expect_equal(midata$n_atoms_index[["3"]], c(2))
})

test_that("atom index is created correctly", {
  atom_index <- create_atom_index(c(1, 1, 3, 4, 4, 3, 6))
  expect_equal(
    atom_index,
    list("1" = c(1, 2), "3" = c(3, 6), "4" = c(4, 5), "6" = c(7))
  )
})

#
# MID normalization
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

test_that("normalizing zero vector gives zero vector", {
  # example "peak area" data, each column an MID
  peak_areas <- matrix(
    c(
      0,
      0,
      0
    ),
    nrow = 3, ncol = 1
  )
  mids <- normalize_mids(peak_areas)
  expect_equal(colSums(mids), c(0))
})

test_that("midata_transform works correctly", {
  midata <- MIData(peak_areas, exp_names)
  # apply identity function
  new_midata <- midata_transform(midata, identity)
  expect_equal(new_midata, midata)
  # apply 13C correction
  new_midata <- midata_transform(midata, c13correct)
  # check corrected MID
  expect_equal(
    get_mids(new_midata, 1, 1)[, 1],
    c13correct(get_mids(midata, 1, 1)[, 1])
  )
  # check averaged corrected MID
  expect_equal(
    get_avg_mid(new_midata, 1, 1),
    c13correct(get_avg_mid(midata, 1, 1))
  )
})
