test_that("simweight gives the same result as the steps one by one", {
  x <- simweight(mouse_pbs_srbc_mini)
  steps <- weight_counts(prepare_counts(mouse_pbs_srbc_mini), threshold = 0.8)
  expect_identical(x, steps)
})

test_that("simweight returns counts, sequences and weighted", {
  x <- simweight(mouse_pbs_srbc_mini)
  expect_true(all(c("counts", "sequences", "weighted") %in% names(x)))
  expect_identical(dim(x$weighted), dim(x$counts))
})

test_that("simweight passes threshold on to weight_counts", {
  x <- simweight(mouse_pbs_srbc_mini, threshold = 0.5)
  steps <- weight_counts(prepare_counts(mouse_pbs_srbc_mini), threshold = 0.5)
  expect_identical(x, steps)
  # and it changes the result compared to the default 0.8
  expect_false(isTRUE(all.equal(x$weighted, simweight(mouse_pbs_srbc_mini)$weighted)))
})

test_that("simweight passes column names on to prepare_counts", {
  d <- mouse_pbs_srbc_mini
  names(d)[names(d) == "cdr3_aa"] <- "junction_aa"
  names(d)[names(d) == "consensus_count"] <- "duplicate_count"
  x <- simweight(d, cdr3 = "junction_aa", count = "duplicate_count")
  expect_identical(x, simweight(mouse_pbs_srbc_mini))
})

test_that("simweight does not depend on nthread", {
  expect_identical(simweight(mouse_pbs_srbc_mini, nthread = 2),
                   simweight(mouse_pbs_srbc_mini))
})

test_that("simweight does not depend on n_cores", {
  # on Windows n_cores > 1 falls back to 1 with a warning
  skip_on_os("windows")
  expect_identical(simweight(mouse_pbs_srbc_mini, n_cores = 2),
                   simweight(mouse_pbs_srbc_mini))
  expect_identical(simweight(mouse_pbs_srbc_mini, method = "blosum", n_cores = 2),
                   simweight(mouse_pbs_srbc_mini, method = "blosum"))
})

test_that("simweight does not filter (yet)", {
  # filtering is left to filter_counts(): every sequence is kept
  x <- simweight(mouse_pbs_srbc_mini)
  n_seqs <- length(unique(paste(mouse_pbs_srbc_mini$v_call, mouse_pbs_srbc_mini$cdr3_aa,
                                mouse_pbs_srbc_mini$j_call)))
  expect_equal(nrow(x$counts), n_seqs)
})
