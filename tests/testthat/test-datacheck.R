test_that("datacheck passes with all required columns", {
  df <- data.frame(
    cdr3_aa = "CASSF",
    consensus_count = 5L,
    v_call = "TRBV1",
    j_call = "TRBJ1",
    sample_processing_id = "S1",
    sequence_id = "seq1",
    cdr3 = "TGTGCC",
    stringsAsFactors = FALSE
  )
  # datacheck always prints status messages; just verify it doesn't error
  expect_no_error(datacheck(df))
})

test_that("datacheck stops when a required column is missing", {
  df <- data.frame(
    cdr3_aa = "CASSF",
    consensus_count = 5L,
    v_call = "TRBV1",
    # j_call missing
    sample_processing_id = "S1",
    sequence_id = "seq1",
    stringsAsFactors = FALSE
  )
  expect_error(datacheck(df), "Dataset incomplete")
})

test_that("datacheck stops when multiple required columns are missing", {
  df <- data.frame(x = 1:3)
  expect_error(datacheck(df), "Dataset incomplete")
})

test_that("datacheck prints a note when cdr3 nucleotide column is absent", {
  df <- data.frame(
    cdr3_aa = "CASSF",
    consensus_count = 5L,
    v_call = "TRBV1",
    j_call = "TRBJ1",
    sample_processing_id = "S1",
    sequence_id = "seq1",
    stringsAsFactors = FALSE
  )
  # Without cdr3 column a note should be printed but no error thrown
  expect_output(datacheck(df), "cdr3 nucleotide sequence not present")
})
