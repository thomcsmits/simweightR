test_that("as_counts_matrix returns a matrix", {
  result <- adjust_counts(mouse_PBSvTCZ_data_minisubset)
  m      <- as_counts_matrix(result)
  expect_true(is.matrix(m))
})

test_that("matrix columns correspond to sample_processing_id values", {
  result  <- adjust_counts(mouse_PBSvTCZ_data_minisubset)
  m       <- as_counts_matrix(result)
  samples <- unique(mouse_PBSvTCZ_data_minisubset$sample_processing_id)
  expect_equal(sort(colnames(m)), sort(samples))
})

test_that("matrix rows correspond to sequence_id values", {
  result <- adjust_counts(mouse_PBSvTCZ_data_minisubset)
  m      <- as_counts_matrix(result)
  seqs   <- unique(mouse_PBSvTCZ_data_minisubset$sequence_id)
  expect_equal(sort(rownames(m)), sort(seqs))
})

test_that("matrix contains no NA values", {
  result <- adjust_counts(mouse_PBSvTCZ_data_minisubset)
  m      <- as_counts_matrix(result)
  expect_false(any(is.na(m)))
})

test_that("matrix values are non-negative", {
  result <- adjust_counts(mouse_PBSvTCZ_data_minisubset)
  m      <- as_counts_matrix(result)
  expect_true(all(m >= 0))
})

test_that("doFilter removes rows present in only one sample", {
  result    <- adjust_counts(mouse_PBSvTCZ_data_minisubset)
  m_all     <- as_counts_matrix(result, doFilter = FALSE)
  m_filtered <- as_counts_matrix(result, doFilter = TRUE)
  # Filtered matrix should have fewer or equal rows than unfiltered
  expect_lte(nrow(m_filtered), nrow(m_all))
  # Every row in filtered matrix must have counts > 0 in at least 2 columns
  if (nrow(m_filtered) > 0) {
    expect_true(all(rowSums(m_filtered > 0) >= 2))
  }
})
