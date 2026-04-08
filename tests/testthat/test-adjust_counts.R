test_that("adjust_counts returns a dataframe with a wrc column", {
  result <- adjust_counts(mouse_PBSvTCZ_data_minisubset)
  expect_true(is.data.frame(result))
  expect_true("wrc" %in% colnames(result))
})

test_that("adjust_counts output has one row per sequence x sample combination", {
  result <- adjust_counts(mouse_PBSvTCZ_data_minisubset)
  # dataprep uses aggregate(drop=FALSE) which creates the full cross-product of
  # all sequence_ids x all sample_processing_ids
  n_seqs <- length(unique(mouse_PBSvTCZ_data_minisubset$sequence_id))
  n_samples <- length(unique(mouse_PBSvTCZ_data_minisubset$sample_processing_id))
  expect_equal(nrow(result), n_seqs * n_samples)
})

test_that("adjust_counts wrc values are non-negative", {
  result <- adjust_counts(mouse_PBSvTCZ_data_minisubset)
  expect_true(all(result$wrc >= 0, na.rm = TRUE))
})

test_that("adjust_counts works with BLOSUM method", {
  result <- adjust_counts(mouse_PBSvTCZ_data_minisubset,
                                   sim_method = "BLOSUM")
  expect_true("wrc" %in% colnames(result))
  expect_true(all(result$wrc >= 0, na.rm = TRUE))
})

test_that("adjust_counts with stricter cutoff does not raise errors", {
  result <- adjust_counts(mouse_PBSvTCZ_data_minisubset, cutoff = 0.99)
  expect_true("wrc" %in% colnames(result))
})

test_that("adjust_counts with permissive cutoff does not raise errors", {
  result <- adjust_counts(mouse_PBSvTCZ_data_minisubset, cutoff = 0.5)
  expect_true("wrc" %in% colnames(result))
})

test_that("adjust_counts errors on input missing required columns", {
  bad_df <- data.frame(x = 1:5)
  expect_error(adjust_counts(bad_df))
})

test_that("adjust_counts accepts a custom similarity function", {
  identity_sim <- function(seq_df) diag(nrow(seq_df))
  result <- adjust_counts(mouse_PBSvTCZ_data_minisubset, sim_method = identity_sim)
  expect_true(is.data.frame(result))
  expect_true("wrc" %in% colnames(result))
})

test_that("custom sim function with identity matrix leaves counts unchanged", {
  # Identity matrix means no neighbours, so wrc should equal consensus_count
  identity_sim <- function(seq_df) diag(nrow(seq_df))
  result <- adjust_counts(mouse_PBSvTCZ_data_minisubset, sim_method = identity_sim)
  expect_equal(result$wrc, result$consensus_count)
})

test_that("adjust_counts errors when custom function returns a non-matrix", {
  bad_sim <- function(seq_df) as.data.frame(diag(nrow(seq_df)))
  expect_error(adjust_counts(mouse_PBSvTCZ_data_minisubset, sim_method = bad_sim),
               "must return a matrix")
})

test_that("adjust_counts errors when custom function returns wrong dimensions", {
  wrong_size_sim <- function(seq_df) diag(nrow(seq_df) + 1)
  expect_error(adjust_counts(mouse_PBSvTCZ_data_minisubset, sim_method = wrong_size_sim),
               regexp = "x")
})

test_that("adjust_counts errors when custom function returns values outside [0, 1]", {
  out_of_range_sim <- function(seq_df) {
    m <- diag(nrow(seq_df))
    m[1, 1] <- 2
    m
  }
  expect_error(adjust_counts(mouse_PBSvTCZ_data_minisubset, sim_method = out_of_range_sim),
               "between 0 and 1")
})

test_that("adjust_counts errors on unrecognised string sim_method", {
  expect_error(adjust_counts(mouse_PBSvTCZ_data_minisubset, sim_method = "UNKNOWN"),
               "HAMMING")
})
