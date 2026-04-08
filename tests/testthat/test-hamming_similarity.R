make_vj_data <- function(seqs, len = nchar(seqs[1])) {
  data.frame(
    cdr3_aa = seqs,
    length = nchar(seqs),
    v_call = "TRBV1",
    j_call = "TRBJ1",
    stringsAsFactors = FALSE
  )
}

test_that("hamming_similarity returns a square matrix", {
  vj <- make_vj_data(c("CASS", "CASG", "CASW"))
  m <- hamming_similarity(vj)
  expect_true(is.matrix(m))
  expect_equal(dim(m), c(3, 3))
})

test_that("hamming_similarity diagonal is 1 (identical sequence)", {
  vj <- make_vj_data(c("CASS", "CASG", "CASW"))
  m <- hamming_similarity(vj)
  expect_equal(diag(m), c(1, 1, 1))
})

test_that("hamming_similarity values are in [0, 1]", {
  vj <- make_vj_data(c("CASS", "CASG", "CASW"))
  m <- hamming_similarity(vj)
  expect_true(all(m >= 0 & m <= 1))
})

test_that("hamming_similarity is symmetric", {
  vj <- make_vj_data(c("CASS", "CASG", "CASW"))
  m <- hamming_similarity(vj)
  # Compare values only; rownames are set but colnames are not, so t(m) has
  # swapped dimnames — strip them before comparing
  expect_equal(unname(m), unname(t(m)))
})

test_that("hamming_similarity computes correct values for known sequences", {
  # "CASS" vs "CASG": 1 mismatch out of 4 -> similarity = 0.75
  vj <- make_vj_data(c("CASS", "CASG"))
  m <- hamming_similarity(vj)
  expect_equal(as.numeric(m[1, 2]), 0.75)
})

test_that("hamming_similarity gives 0 for completely different sequences", {
  # All 4 positions differ
  vj <- make_vj_data(c("CASS", "WQLP"))
  m <- hamming_similarity(vj)
  expect_equal(as.numeric(m[1, 2]), 0)
})

test_that("hamming_similarity errors with fewer than 2 sequences", {
  vj <- make_vj_data(c("CASS"))
  expect_error(hamming_similarity(vj), "At least two sequences")
})

test_that("hamming_similarity errors when sequences have different lengths", {
  vj <- data.frame(
    cdr3_aa = c("CASS", "CASSG"),
    length = c(4, 5),
    v_call = "TRBV1",
    j_call = "TRBJ1",
    stringsAsFactors = FALSE
  )
  expect_error(hamming_similarity(vj), "same length")
})
