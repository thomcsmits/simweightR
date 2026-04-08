test_that("adjust_counts returns a dataframe with a wrc column", {
  result <- adjust_counts(mouse_PBSvTCZ_data_minisubset)
  expect_true(is.data.frame(result))
  expect_true("wrc" %in% colnames(result))
})

test_that("adjust_counts output has one row per sequence x sample combination", {
  result <- adjust_counts(mouse_PBSvTCZ_data_minisubset)
  # dataprep uses aggregate(drop=FALSE) which creates the full cross-product of
  # all sequence_ids x all sample_processing_ids
  n_seqs    <- length(unique(mouse_PBSvTCZ_data_minisubset$sequence_id))
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
