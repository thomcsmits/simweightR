test_that("update_counts returns same length vector as input", {
  counts <- c(10, 5, 3)
  sim    <- diag(3)
  result <- update_counts(counts, sim, cutoff = 0.8)
  expect_equal(length(result), 3)
})

test_that("update_counts leaves counts unchanged when no neighbors exceed cutoff", {
  counts <- c(10, 5, 3)
  # Identity matrix: each sequence is only similar to itself (sim = 1 >= cutoff)
  sim <- diag(3)
  result <- update_counts(counts, sim, cutoff = 0.8)
  # For count[i] > 0: weighted_sum = count[i] + sim[i,i]*count[i] = 2*count[i]
  #                   total_weight = 1 + sim[i,i] = 2
  #                   result = count[i]  (self-similarity cancels out)
  expect_equal(result, counts)
})

test_that("update_counts incorporates neighbor counts above cutoff", {
  # Two identical sequences (sim = 1); both influence each other.
  # Note: similar_indices includes the sequence itself (sim[i,i]=1 >= cutoff),
  # so the self-count is double-weighted in the sum.
  counts <- c(10, 0)
  sim    <- matrix(c(1, 1, 1, 1), nrow = 2)
  result <- update_counts(counts, sim, cutoff = 0.8)

  # seq1 (count=10 > 0): similar_indices = {1, 2}
  #   weighted_sum  = 10 + sim[1,1]*10 + sim[1,2]*0 = 20
  #   total_weight  = 1 + sim[1,1] + sim[1,2] = 3
  #   result[1] = 20/3
  # seq2 (count=0): similar_indices = {1, 2}
  #   weighted_sum  = sim[2,1]*10 + sim[2,2]*0 = 10
  #   total_weight  = sim[2,1] + sim[2,2] = 2
  #   result[2] = 10/2 = 5
  expect_equal(result[1], 20 / 3)
  expect_equal(result[2], 5)
})

test_that("update_counts returns 0 for zero-count sequence with no similar neighbors", {
  counts <- c(5, 0)
  sim    <- matrix(c(1, 0, 0, 1), nrow = 2)  # no cross-similarity
  result <- update_counts(counts, sim, cutoff = 0.8)
  expect_equal(result[2], 0)
})

test_that("update_counts respects the cutoff threshold", {
  counts <- c(10, 8)
  # sim[1,2] = 0.7, below cutoff of 0.8 -> sequences do not influence each other
  sim    <- matrix(c(1, 0.7, 0.7, 1), nrow = 2)
  result <- update_counts(counts, sim, cutoff = 0.8)
  expect_equal(result, counts)
})

test_that("update_counts weighted average is correct for two equal sequences", {
  counts <- c(8, 4)
  # Perfect similarity: all sim = 1; self is included in similar_indices
  sim    <- matrix(c(1, 1, 1, 1), nrow = 2)
  result <- update_counts(counts, sim, cutoff = 0.8)

  # seq1: ws = 8 + 1*8 + 1*4 = 20; tw = 1+1+1 = 3; result = 20/3
  # seq2: ws = 4 + 1*4 + 1*8 = 16; tw = 3;          result = 16/3
  expect_equal(result[1], 20 / 3)
  expect_equal(result[2], 16 / 3)
})
