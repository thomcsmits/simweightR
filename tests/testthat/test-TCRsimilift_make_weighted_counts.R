test_that("TCRsimilift_make_weighted_counts returns a matrix", {
  result <- TCRsimilift_calculate(mouse_PBSvTCZ_data_minisubset)
  m      <- TCRsimilift_make_weighted_counts(result)
  expect_true(is.matrix(m))
})

test_that("matrix columns correspond to sample_processing_id values", {
  result  <- TCRsimilift_calculate(mouse_PBSvTCZ_data_minisubset)
  m       <- TCRsimilift_make_weighted_counts(result)
  samples <- unique(mouse_PBSvTCZ_data_minisubset$sample_processing_id)
  expect_equal(sort(colnames(m)), sort(samples))
})

test_that("matrix rows correspond to sequence_id values", {
  result <- TCRsimilift_calculate(mouse_PBSvTCZ_data_minisubset)
  m      <- TCRsimilift_make_weighted_counts(result)
  seqs   <- unique(mouse_PBSvTCZ_data_minisubset$sequence_id)
  expect_equal(sort(rownames(m)), sort(seqs))
})

test_that("matrix contains no NA values", {
  result <- TCRsimilift_calculate(mouse_PBSvTCZ_data_minisubset)
  m      <- TCRsimilift_make_weighted_counts(result)
  expect_false(any(is.na(m)))
})

test_that("matrix values are non-negative", {
  result <- TCRsimilift_calculate(mouse_PBSvTCZ_data_minisubset)
  m      <- TCRsimilift_make_weighted_counts(result)
  expect_true(all(m >= 0))
})

test_that("doFilter removes rows present in only one sample", {
  result    <- TCRsimilift_calculate(mouse_PBSvTCZ_data_minisubset)
  m_all     <- TCRsimilift_make_weighted_counts(result, doFilter = FALSE)
  m_filtered <- TCRsimilift_make_weighted_counts(result, doFilter = TRUE)
  # Filtered matrix should have fewer or equal rows than unfiltered
  expect_lte(nrow(m_filtered), nrow(m_all))
  # Every row in filtered matrix must have counts > 0 in at least 2 columns
  if (nrow(m_filtered) > 0) {
    expect_true(all(rowSums(m_filtered > 0) >= 2))
  }
})
