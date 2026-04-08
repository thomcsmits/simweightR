test_that("TCRsimilift_calculate returns a dataframe with a wrc column", {
  result <- TCRsimilift_calculate(mouse_PBSvTCZ_data_minisubset)
  expect_true(is.data.frame(result))
  expect_true("wrc" %in% colnames(result))
})

test_that("TCRsimilift_calculate output has one row per sequence x sample combination", {
  result <- TCRsimilift_calculate(mouse_PBSvTCZ_data_minisubset)
  # dataprep uses aggregate(drop=FALSE) which creates the full cross-product of
  # all sequence_ids x all sample_processing_ids
  n_seqs    <- length(unique(mouse_PBSvTCZ_data_minisubset$sequence_id))
  n_samples <- length(unique(mouse_PBSvTCZ_data_minisubset$sample_processing_id))
  expect_equal(nrow(result), n_seqs * n_samples)
})

test_that("TCRsimilift_calculate wrc values are non-negative", {
  result <- TCRsimilift_calculate(mouse_PBSvTCZ_data_minisubset)
  expect_true(all(result$wrc >= 0, na.rm = TRUE))
})

test_that("TCRsimilift_calculate works with BLOSUM method", {
  result <- TCRsimilift_calculate(mouse_PBSvTCZ_data_minisubset,
                                   sim_method = "BLOSUM")
  expect_true("wrc" %in% colnames(result))
  expect_true(all(result$wrc >= 0, na.rm = TRUE))
})

test_that("TCRsimilift_calculate with stricter cutoff does not raise errors", {
  result <- TCRsimilift_calculate(mouse_PBSvTCZ_data_minisubset, cutoff = 0.99)
  expect_true("wrc" %in% colnames(result))
})

test_that("TCRsimilift_calculate with permissive cutoff does not raise errors", {
  result <- TCRsimilift_calculate(mouse_PBSvTCZ_data_minisubset, cutoff = 0.5)
  expect_true("wrc" %in% colnames(result))
})

test_that("TCRsimilift_calculate errors on input missing required columns", {
  bad_df <- data.frame(x = 1:5)
  expect_error(TCRsimilift_calculate(bad_df))
})
