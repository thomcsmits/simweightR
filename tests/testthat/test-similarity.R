test_that("hamming_similarity takes a character vector", {
  vj <- data.frame(cdr3_aa = (c("CASS", "CASG", "CASW")))
  expect_equal(hamming_similarity(vj$cdr3_aa), hamming_similarity(c("CASS", "CASG", "CASW")))
})

test_that("hamming_similarity errors when input is not a character vector", {
  vj <- data.frame(cdr3_aa = c("CASS", "CASG"))
  expect_error(hamming_similarity(vj), "character vector")
  expect_error(hamming_similarity(vj["cdr3_aa"]), "character vector")
})

test_that("hamming_similarity errors when sequences have different lengths", {
  expect_error(hamming_similarity(c("CASS", "CASSG")), "same length")
})

test_that("hamming_similarity returns a square matrix", {
  m <- hamming_similarity(c("CASS", "CASG", "CASW"))
  expect_true(is.matrix(m))
  expect_equal(dim(m), c(3, 3))
})

test_that("hamming_similarity diagonal is 1 (identical sequence)", {
  m <- hamming_similarity(c("CASS", "CASG", "CASW"))
  expect_equal(diag(m), c(1, 1, 1))
})

test_that("hamming_similarity values are in [0, 1]", {
  m <- hamming_similarity(c("CASS", "CASG", "CASW"))
  expect_true(all(m >= 0 & m <= 1))
})

test_that("hamming_similarity is symmetric", {
  m <- hamming_similarity(c("CASS", "CASG", "CASW"))
  expect_equal(m, t(m))
})

test_that("hamming_similarity computes correct values for known sequences", {
  # "CASS" vs "CASG": 1 mismatch out of 4 -> similarity = 0.75
  m <- hamming_similarity(c("CASS", "CASG"))
  expect_equal(m[1, 2], 0.75)
})

test_that("hamming_similarity gives 0 for completely different sequences", {
  # All 4 positions differ
  m <- hamming_similarity(c("CASS", "WQLP"))
  expect_equal(m[1, 2], 0)
})

# reference: count mismatching positions directly
hamming_reference <- function(sequences) {
  n <- length(sequences)
  chars <- strsplit(sequences, "")
  sim <- matrix(0, n, n)
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      sim[i, j] <- 1 - sum(chars[[i]] != chars[[j]]) / length(chars[[i]])
    }
  }
  sim
}

test_that("hamming_similarity matches direct count", {
  s <- c("CASSLGQETQYF", "CASSLGRETQYF", "CASSPGQETQYA", "CASRLGQGTQYF")
  expect_equal(hamming_similarity(s), hamming_reference(s))
})

test_that("hamming_similarity gives 1 for duplicate sequences", {
  s <- c("CASSLGQETQYF", "CASSLGRETQYF", "CASSLGQETQYF")
  m <- hamming_similarity(s)
  expect_equal(m[1, 3], 1)
  expect_equal(m[1, ], m[3, ])
})

test_that("hamming_similarity handles a single sequence", {
  m <- hamming_similarity("CASS")
  expect_equal(m, matrix(1))
})

test_that("blosum_similarity computes correct values for known sequences", {
  # "CAS" vs "CAT": S -> T at position 3, d(S, T) = 0.0985
  #   similarity = 1 - 0.0985 / 3 = 0.9672
  # "CAS" vs "CIS": A -> I at position 2, d(A, I) = 0.2978
  #   similarity = 1 - 0.2978 / 3 = 0.9007
  # "CAT" vs "CIS": both differences, 1 - (0.0985 + 0.2978) / 3 = 0.8679
  m <- blosum_similarity(c("CAS", "CAT", "CIS"))
  expect_equal(m[1, 2], 1 - 0.0985 / 3)
  expect_equal(m[1, 3], 1 - 0.2978 / 3)
  expect_equal(m[2, 3], 1 - (0.0985 + 0.2978) / 3)
  # a conservative substitution (S -> T) is more similar than A -> I
  expect_gt(m[1, 2], m[1, 3])
})

test_that("blosum_similarity returns a symmetric matrix with 1 on the diagonal", {
  # includes a duplicate sequence (1 and 4)
  s <- c("CASSLGQETQYF", "CASSLGRETQYF", "CASSIGQETQYF", "CASSLGQETQYF")
  m <- blosum_similarity(s)
  expect_true(is.matrix(m))
  expect_equal(dim(m), c(4, 4))
  expect_equal(m, t(m))
  expect_equal(diag(m), rep(1, 4))
  expect_true(all(m >= 0 & m <= 1))
  # duplicates have similarity 1 and the same row
  expect_equal(m[1, 4], 1)
  expect_equal(m[1, ], m[4, ])
})

# reference: for each pair, look up the distance at each position by name
blosum_reference <- function(sequences) {
  n <- length(sequences)
  chars <- strsplit(sequences, "")
  sim <- matrix(0, n, n)
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      d <- simweightR::blosum62_distance[cbind(chars[[i]], chars[[j]])]
      sim[i, j] <- 1 - sum(d) / length(chars[[i]])
    }
  }
  sim
}

test_that("blosum_similarity matches per-pair lookup", {
  x <- prepare_counts(mouse_pbs_srbc_mini)
  s <- x$sequences$cdr3_aa[x$sequences$j_call == "TRBJ1-1"]
  expect_equal(blosum_similarity(s), blosum_reference(s))
})

test_that("blosum_similarity with a 0/1 distance equals hamming_similarity", {
  # Hamming is the special case where the distance is 0 for equal and 1 for
  # different amino acids
  zero_one <- 1 - diag(20)
  dimnames(zero_one) <- dimnames(blosum62_distance)
  s <- c("CASSLGQETQYF", "CASSLGRETQYF", "CASSPGQETQYA", "CASRLGQGTQYF")
  expect_equal(blosum_similarity(s, distance = zero_one), hamming_similarity(s))
})

test_that("blosum_similarity handles a single sequence", {
  expect_equal(blosum_similarity("CASS"), matrix(1))
})

test_that("blosum_similarity rejects bad input", {
  expect_error(blosum_similarity(c("CASS", "CASSG")), "same length")
  vj <- data.frame(cdr3_aa = c("CASS", "CASG"))
  expect_error(blosum_similarity(vj), "character vector")
  # X is not one of the 20 amino acids in blosum62_distance
  expect_error(blosum_similarity(c("CASX", "CASS")), "not in the distance matrix")
})

test_that("check_similarity accepts valid similarity matrices", {
  s <- c("CASSLGQETQYF", "CASSLGRETQYF", "CASSPGQETQYA")
  expect_no_error(check_similarity(hamming_similarity(s), 3))
  expect_no_error(check_similarity(blosum_similarity(s), 3))
  # a single sequence
  expect_no_error(check_similarity(matrix(1), 1))
  # values exactly 0 and 1 are allowed
  expect_no_error(check_similarity(matrix(c(1, 0, 0, 1), 2), 2))
})

test_that("check_similarity rejects invalid similarity matrices", {
  # not a matrix
  expect_error(check_similarity(as.data.frame(diag(2)), 2), "must return a matrix")
  expect_error(check_similarity(c(1, 1), 2), "must return a matrix")
  # wrong dimensions, also when not square
  expect_error(check_similarity(diag(3), 2), "must be 2 x 2, got 3 x 3")
  expect_error(check_similarity(matrix(1, 2, 3), 2), "must be 2 x 2, got 2 x 3")
  # values outside [0, 1]
  expect_error(check_similarity(matrix(c(1, 2, 2, 1), 2), 2), "between 0 and 1")
  expect_error(check_similarity(matrix(c(1, -0.1, -0.1, 1), 2), 2), "between 0 and 1")
})
