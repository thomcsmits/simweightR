test_that("weighting returns matrix with same dimensions as counts", {
  counts <- cbind(c(10, 5, 3), c(1, 2, 3))
  sim <- diag(3)
  result <- weight_group(counts, sim, threshold = 0.8)
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(3, 2))
})

test_that("weighting leaves counts unchanged when there are no neighbours", {
  counts <- cbind(c(10, 5, 3))
  # Identity matrix: each sequence is only similar to itself
  sim <- diag(3)
  # For each i: weighted_sum = sim[i,i]*count[i] = count[i]
  #             total_weight = sim[i,i] = 1
  #             result = count[i]  (self is the only neighbour)
  result <- weight_group(counts, sim, threshold = 0.8)
  expect_equal(result, counts)
})

test_that("weighting incorporates neighbour counts above threshold", {
  # Two identical sequences (sim = 1); both influence each other.
  # Note: the sequence itself is included once (sim[i,i] = 1), for every
  # sequence, whether its count is 0 or not.
  counts <- cbind(c(10, 0))
  sim <- matrix(c(1, 1, 1, 1), nrow = 2)
  result <- weight_group(counts, sim, threshold = 0.8)

  # seq1 (count=10): neighbours = {1, 2}
  #   weighted_sum  = sim[1,1]*10 + sim[1,2]*0 = 10
  #   total_weight  = sim[1,1] + sim[1,2] = 2
  #   result[1] = 10/2 = 5
  # seq2 (count=0): neighbours = {1, 2}
  #   weighted_sum  = sim[2,1]*10 + sim[2,2]*0 = 10
  #   total_weight  = sim[2,1] + sim[2,2] = 2
  #   result[2] = 10/2 = 5
  expect_equal(result[1, 1], 5)
  expect_equal(result[2, 1], 5)
})

test_that("weighting returns 0 for zero-count sequence with no similar neighbours", {
  counts <- cbind(c(5, 0))
  sim <- matrix(c(1, 0, 0, 1), nrow = 2)  # no cross-similarity
  result <- weight_group(counts, sim, threshold = 0.8)
  expect_equal(result[2, 1], 0)
})

test_that("weighting respects the threshold", {
  counts <- cbind(c(10, 8))
  # sim[1,2] = 0.7, below threshold of 0.8 -> sequences do not influence each other
  sim <- matrix(c(1, 0.7, 0.7, 1), nrow = 2)
  result <- weight_group(counts, sim, threshold = 0.8)
  expect_equal(result, counts)
})

test_that("weighting weighted average is correct for two equal sequences", {
  counts <- cbind(c(8, 4))
  # Perfect similarity: all sim = 1; self is included once
  sim <- matrix(c(1, 1, 1, 1), nrow = 2)
  result <- weight_group(counts, sim, threshold = 0.8)

  # seq1: ws = 1*8 + 1*4 = 12; tw = 1+1 = 2; result = 12/2 = 6
  # seq2: ws = 1*4 + 1*8 = 12; tw = 2;       result = 12/2 = 6
  expect_equal(result[1, 1], 6)
  expect_equal(result[2, 1], 6)
})

test_that("weighting gives the toy example values", {
  # seq1 and seq2 are neighbours (6/7 >= 0.8); seq3 is not (0.5 < 0.8)
  counts <- cbind(c(10, 0, 4))
  sim <- matrix(c(1, 6 / 7, 0.5,
                  6 / 7, 1, 0.5,
                  0.5, 0.5, 1), nrow = 3)
  result <- weight_group(counts, sim, threshold = 0.8)

  # seq1: ws = 1*10 + 6/7*0 = 10;  tw = 1 + 6/7 = 13/7; result = 70/13 = 5.38
  # seq2: ws = 6/7*10 + 1*0 = 60/7; tw = 13/7;          result = 60/13 = 4.62
  # seq3: ws = 1*4 = 4;             tw = 1;             result = 4
  expect_equal(result[, 1], c(70 / 13, 60 / 13, 4))
})

test_that("weighting weights each sample independently", {
  # weighting all samples at once gives the same as one sample at a time
  counts <- cbind(c(10, 0, 4), c(1, 2, 3), c(0, 0, 7))
  sim <- matrix(c(1, 6 / 7, 0.5,
                  6 / 7, 1, 0.5,
                  0.5, 0.5, 1), nrow = 3)
  result <- weight_group(counts, sim, threshold = 0.8)
  for (s in seq_len(ncol(counts))) {
    expect_equal(result[, s], weight_group(counts[, s, drop = FALSE], sim, threshold = 0.8)[, 1])
  }
})

test_that("weighting counts pairs exactly at the threshold as neighbours", {
  # 1 - 3/20 is 0.85 mathematically, but stored as 0.84999999999999997780,
  # while seq(0.5, 0.95, by = 0.05)[8] is stored as 0.85000000000000008882
  s <- 1 - 3 / 20
  threshold <- seq(0.5, 0.95, by = 0.05)[8]
  expect_false(s >= threshold)  # plain comparison fails

  counts <- cbind(c(10, 0))
  sim <- matrix(c(1, s, s, 1), nrow = 2)
  result <- weight_group(counts, sim, threshold = threshold)

  # seq1: ws = 1*10 + s*0 = 10; tw = 1 + s; result = 10 / (1 + s)
  # seq2: ws = s*10 + 1*0 = 10s; tw = 1 + s; result = 10s / (1 + s)
  expect_equal(result[, 1], c(10 / (1 + s), 10 * s / (1 + s)))
})

test_that("weighting warns and sets the diagonal to 1 if it is not 1", {
  counts <- cbind(c(10, 0, 4))
  sim <- matrix(c(1, 6 / 7, 0.5,
                  6 / 7, 1, 0.5,
                  0.5, 0.5, 1), nrow = 3)
  sim0 <- sim
  diag(sim0) <- 0
  # with diagonal 0 the result is the same as with diagonal 1 (toy example)
  expect_warning(result <- weight_group(counts, sim0, threshold = 0.8), "diagonal")
  expect_equal(result, weight_group(counts, sim, threshold = 0.8))
})

test_that("weighting rejects bad input", {
  sim <- diag(3)
  # counts not a matrix
  expect_error(weight_group(c(10, 0, 4), sim, threshold = 0.8), "must be a matrix")
  # number of rows does not match
  expect_error(weight_group(cbind(c(10, 0)), sim, threshold = 0.8), "one row per sequence")
  # threshold outside [0, 1]
  expect_error(weight_group(cbind(c(10, 0, 4)), sim, threshold = 1.5), "between 0 and 1")
  expect_error(weight_group(cbind(c(10, 0, 4)), sim, threshold = -0.1), "between 0 and 1")
})
