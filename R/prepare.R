#' AIRR table to counts matrix & sequences info
#'
#' Get counts per sample for V gene + J gene + CDR3 amino acid sequence group
#'
#' @inheritParams simweight
#'
#' @returns List with
#'   * sequences: data frame with id, cdr3_aa, v_call and j_call, one row per
#'     sequence
#'
#'   * counts: matrix of counts, one row per sequence (rownames "V_CDR3_J") and
#'     one column per sample, in the same order as the rows of sequences
#'
#' @export
prepare_counts <- function(df,
                           cdr3 = "cdr3_aa",
                           v = "v_call",
                           j = "j_call",
                           sample = "sample_processing_id",
                           count = "consensus_count") {
  required <- c(cdr3, v, j, sample, count)
  missing <- setdiff(required, colnames(df))
  if (length(missing) > 0) {
    stop("df is missing columns: ", paste(missing, collapse = ", "))
  }

  id <- paste(df[[v]], df[[cdr3]], df[[j]], sep = "_")
  # get levels for factors, order of first appearance, so counts sorting does
  # not depend on locale-specific sorting
  ids <- unique(id)
  samples <- unique(df[[sample]])

  # sum counts per sequence and sample (missing = 0)
  # matrix cell of each row: (row + (column - 1) * n_rows)
  cell <- (match(df[[sample]], samples) - 1) * length(ids) + match(id, ids)
  # sum counts per cell (rows of df that are the same sequence in the same
  # sample)
  sums <- rowsum(df[[count]], cell, reorder = FALSE)
  # other cells (sequence not observed in that sample) stay 0
  counts <- matrix(0, length(ids), length(samples),
                   dimnames = list(ids, samples))
  counts[as.numeric(rownames(sums))] <- sums

  # fix column names
  first <- match(ids, id)
  sequences <- data.frame(id = ids,
                          cdr3_aa = df[[cdr3]][first],
                          v_call = df[[v]][first],
                          j_call = df[[j]][first])

  list(sequences = sequences, counts = counts)
}
