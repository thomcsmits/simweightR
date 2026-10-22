# toy AIRR table: rows 1, 3 and 5 are the same sequence (V1_CASSA_J1); row 4
# has the same CDR3 but a different V gene, so it is a different sequence
make_toy <- function() {
  data.frame(cdr3_aa = c("CASSA", "CASSB", "CASSA", "CASSA", "CASSA"),
             v_call = c("V1", "V1", "V1", "V2", "V1"),
             j_call = "J1",
             sample_processing_id = c("s2", "s1", "s1", "s1", "s2"),
             consensus_count = c(5, 3, 2, 7, 1))
}

test_that("prepare_counts gives the toy example matrix", {
  x <- prepare_counts(make_toy())

  # V1_CASSA_J1 in s2: rows 1 + 5 = 5 + 1 = 6; in s1: row 3 = 2
  # V1_CASSB_J1 not observed in s2 -> 0
  # rows and columns in order of first appearance (not sorted)
  expected <- matrix(c(6, 0, 0,
                       2, 3, 7), nrow = 3,
                     dimnames = list(c("V1_CASSA_J1", "V1_CASSB_J1", "V2_CASSA_J1"),
                                     c("s2", "s1")))
  expect_equal(x$counts, expected)
})

test_that("prepare_counts keeps the same CDR3 with a different V gene separate", {
  x <- prepare_counts(make_toy())
  expect_true(all(c("V1_CASSA_J1", "V2_CASSA_J1") %in% rownames(x$counts)))
  expect_equal(x$counts["V2_CASSA_J1", "s1"], 7)
})

test_that("prepare_counts sequences table lines up with the rows of counts", {
  x <- prepare_counts(make_toy())
  expect_identical(x$sequences$id, rownames(x$counts))
  expect_equal(x$sequences$cdr3_aa, c("CASSA", "CASSB", "CASSA"))
  expect_equal(x$sequences$v_call, c("V1", "V1", "V2"))
  expect_equal(x$sequences$j_call, c("J1", "J1", "J1"))
})

test_that("prepare_counts preserves the total count", {
  x <- prepare_counts(mouse_pbs_srbc_mini)
  expect_equal(sum(x$counts), sum(mouse_pbs_srbc_mini$consensus_count))
})

test_that("prepare_counts accepts other column names", {
  # e.g. junction_aa, repertoire_id and duplicate_count instead of the defaults
  d <- make_toy()
  names(d) <- c("junction_aa", "v", "j", "repertoire_id", "duplicate_count")
  x <- prepare_counts(d, cdr3 = "junction_aa", v = "v", j = "j",
                      sample = "repertoire_id", count = "duplicate_count")
  expect_equal(x, prepare_counts(make_toy()))
  # output always uses the fixed names
  expect_named(x$sequences, c("id", "cdr3_aa", "v_call", "j_call"))
})

test_that("prepare_counts errors on input missing required columns", {
  d <- make_toy()
  d$v_call <- NULL
  expect_error(prepare_counts(d), "missing columns: v_call")
  # none of the required columns
  bad_df <- data.frame(x = 1:5)
  expect_error(prepare_counts(bad_df), "missing columns")
})
