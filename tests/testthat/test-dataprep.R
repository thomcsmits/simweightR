make_minimal_df <- function() {
  data.frame(
    sequence_id = c("seq1", "seq1", "seq2"),
    sample_processing_id = c("S1", "S1", "S1"),
    consensus_count = c(3L, 2L, 10L),
    cdr3_aa = c("CASS", "CASS", "CASW"),
    v_call = c("TRBV1", "TRBV1", "TRBV2"),
    j_call = c("TRBJ1", "TRBJ1", "TRBJ2"),
    stringsAsFactors = FALSE
  )
}

test_that("dataprep returns a dataframe with expected columns", {
  result <- dataprep(make_minimal_df())
  expected_cols <- c("sequence_id", "sample_processing_id", "consensus_count",
                     "length", "v_call", "j_call", "cdr3_aa")
  expect_true(all(expected_cols %in% colnames(result)))
})

test_that("dataprep adds correct length column", {
  result <- dataprep(make_minimal_df())
  lengths <- result$length[match(result$sequence_id, result$sequence_id)]
  expect_true(all(result$length == nchar(result$cdr3_aa)))
})

test_that("dataprep aggregates duplicate sequence_id + sample_processing_id counts", {
  result <- dataprep(make_minimal_df())
  # seq1 in S1 has two rows with counts 3 and 2; should sum to 5
  seq1_count <- result$consensus_count[result$sequence_id == "seq1" &
                                         result$sample_processing_id == "S1"]
  expect_equal(seq1_count, 5)
})

test_that("dataprep produces one row per sequence_id x sample_processing_id", {
  df <- make_minimal_df()
  result <- dataprep(df)
  n_unique <- nrow(unique(df[c("sequence_id", "sample_processing_id")]))
  expect_equal(nrow(result), n_unique)
})

test_that("dataprep replaces NA counts with 0", {
  # aggregate with drop=FALSE creates NAs for missing combinations
  df <- data.frame(
    sequence_id = c("seq1", "seq2"),
    sample_processing_id = c("S1", "S2"),
    consensus_count = c(5L, 3L),
    cdr3_aa = c("CASS", "CASW"),
    v_call = c("TRBV1", "TRBV2"),
    j_call = c("TRBJ1", "TRBJ2"),
    stringsAsFactors = FALSE
  )
  result <- dataprep(df)
  expect_false(any(is.na(result$consensus_count)))
})
