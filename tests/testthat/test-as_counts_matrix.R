test_that("as_counts_matrix returns a matrix", {
  result <- adjust_counts(mouse_PBSvTCZ_data_minisubset)
  m <- as_counts_matrix(result)
  expect_true(is.matrix(m))
})

test_that("matrix columns correspond to sample_processing_id values", {
  result <- adjust_counts(mouse_PBSvTCZ_data_minisubset)
  m <- as_counts_matrix(result)
  samples <- unique(mouse_PBSvTCZ_data_minisubset$sample_processing_id)
  expect_equal(sort(colnames(m)), sort(samples))
})

test_that("matrix rows correspond to sequence_id values", {
  result <- adjust_counts(mouse_PBSvTCZ_data_minisubset)
  m <- as_counts_matrix(result)
  seqs <- unique(mouse_PBSvTCZ_data_minisubset$sequence_id)
  expect_equal(sort(rownames(m)), sort(seqs))
})

test_that("matrix contains no NA values", {
  result <- adjust_counts(mouse_PBSvTCZ_data_minisubset)
  m <- as_counts_matrix(result)
  expect_false(any(is.na(m)))
})

test_that("matrix values are non-negative", {
  result <- adjust_counts(mouse_PBSvTCZ_data_minisubset)
  m <- as_counts_matrix(result)
  expect_true(all(m >= 0))
})

test_that("doFilter removes rows present in only one sample", {
  result <- adjust_counts(mouse_PBSvTCZ_data_minisubset)
  m_all <- as_counts_matrix(result, doFilter = FALSE)
  m_filtered <- as_counts_matrix(result, doFilter = TRUE)
  # Filtered matrix should have fewer or equal rows than unfiltered
  expect_lte(nrow(m_filtered), nrow(m_all))
  # Every row in filtered matrix must have counts > 0 in at least 2 columns
  if (nrow(m_filtered) > 0) {
    expect_true(all(rowSums(m_filtered > 0) >= 2))
  }
})

test_that("weighted = FALSE returns consensus_count values", {
  result <- adjust_counts(mouse_PBSvTCZ_data_minisubset)
  m_unweighted <- as_counts_matrix(result, weighted = FALSE)
  # unweighted matrix values should match consensus_count from the result
  for (seq in rownames(m_unweighted)) {
    for (samp in colnames(m_unweighted)) {
      row <- result[result$sequence_id == seq & result$sample_processing_id == samp, ]
      expected <- if (nrow(row) == 0) 0 else row$consensus_count
      expect_equal(m_unweighted[seq, samp], expected)
    }
  }
})

test_that("weighted = FALSE and weighted = TRUE differ when counts were adjusted", {
  result <- adjust_counts(mouse_PBSvTCZ_data_minisubset)
  m_weighted <- as_counts_matrix(result, weighted = TRUE)
  m_unweighted <- as_counts_matrix(result, weighted = FALSE)
  # If any sequences were adjusted, the matrices should differ
  if (any(result$wrc != result$consensus_count)) {
    expect_false(identical(m_weighted, m_unweighted))
  }
})

test_that("weighted = FALSE with doFilter still filters on consensus_count", {
  result <- adjust_counts(mouse_PBSvTCZ_data_minisubset)
  m_filtered <- as_counts_matrix(result, doFilter = TRUE, weighted = FALSE)
  expect_true(is.matrix(m_filtered))
  if (nrow(m_filtered) > 0) {
    expect_true(all(rowSums(m_filtered > 0) >= 2))
  }
})
