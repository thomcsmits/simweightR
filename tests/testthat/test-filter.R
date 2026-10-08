# toy AIRR table:
#   V1_CASSA_J1: s1 = 6, s2 = 2         -> observed in 2 samples
#   V1_CASSB_J1: s2 = 3                 -> observed in 1 sample
#   V2_CASSA_J1: s3 = 40                -> observed in 1 sample (high count)
make_toy <- function() {
  data.frame(cdr3_aa = c("CASSA", "CASSA", "CASSB", "CASSA"),
             v_call = c("V1", "V1", "V1", "V2"),
             j_call = "J1",
             sample_processing_id = c("s1", "s2", "s2", "s3"),
             consensus_count = c(6, 2, 3, 40))
}

test_that("filter_counts keeps sequences observed in at least 2 samples", {
  x <- filter_counts(prepare_counts(make_toy()))
  expect_equal(rownames(x$counts), "V1_CASSA_J1")
  expect_equal(x$sequences$id, "V1_CASSA_J1")
  # columns (samples) are not removed
  expect_equal(colnames(x$counts), c("s1", "s2", "s3"))
})

test_that("filter_counts respects min_samples", {
  x <- prepare_counts(make_toy())
  # min_samples = 1 keeps every sequence
  expect_equal(nrow(filter_counts(x, min_samples = 1)$counts), 3)
  # min_samples = 3: no sequence is observed in all 3 samples
  expect_equal(nrow(filter_counts(x, min_samples = 3)$counts), 0)
})

test_that("filter_counts respects min_count", {
  x <- prepare_counts(make_toy())
  # with min_count = 3, V1_CASSA_J1 (6 and 2) counts as observed in s1 only,
  # so with min_samples = 1 only the sequences with a count >= 3 remain
  kept <- filter_counts(x, min_samples = 1, min_count = 3)
  expect_equal(rownames(kept$counts), c("V1_CASSA_J1", "V1_CASSB_J1", "V2_CASSA_J1"))
  kept <- filter_counts(x, min_samples = 1, min_count = 10)
  expect_equal(rownames(kept$counts), "V2_CASSA_J1")
  # and V1_CASSA_J1 is no longer in 2 samples
  expect_equal(nrow(filter_counts(x, min_samples = 2, min_count = 3)$counts), 0)
})

test_that("filter_counts keeps counts, sequences and weighted aligned", {
  x <- weight_counts(prepare_counts(mouse_pbs_srbc_mini), threshold = 0.8)
  kept <- filter_counts(x)
  expect_identical(kept$sequences$id, rownames(kept$counts))
  expect_identical(rownames(kept$weighted), rownames(kept$counts))
  # the kept rows are unchanged
  expect_identical(kept$weighted, x$weighted[rownames(kept$counts), , drop = FALSE])
  expect_identical(kept$counts, x$counts[rownames(kept$counts), , drop = FALSE])
})

test_that("filter_counts works without weighted", {
  x <- filter_counts(prepare_counts(make_toy()))
  expect_null(x$weighted)
})

test_that("filter_counts is judged on the unweighted counts", {
  # CASSA and CASSB are neighbours (similarity 4/5 = 0.8), each observed in
  # one sample only. After weighting, both have a weighted count > 0 in both
  # samples, but they are still removed.
  d <- data.frame(cdr3_aa = c("CASSA", "CASSB"), v_call = "V1", j_call = "J1",
                  sample_processing_id = c("s1", "s2"), consensus_count = c(10, 5))
  x <- weight_counts(prepare_counts(d), threshold = 0.8)
  expect_true(all(x$weighted > 0))
  expect_equal(nrow(filter_counts(x)$counts), 0)
})
