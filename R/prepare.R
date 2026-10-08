#' AIRR table to counts matrix & sequences info
#'
#' Get counts per sample for V gene + J gene + CD3 amino acid sequence group
#' 
#' @param df AIRR formatted data frame with one row per sequence
#' @param cdr3 Column with the CDR3 amino acid sequence, e.g. "cdr3_aa" or "junction_aa" (junction = CDR3 plus the conserved flanking C and F/W)
#' @param v Column with the V gene
#' @param j Column with the J gene
#' @param sample Column with the sample identifier, e.g.
#'   "sample_processing_id" or "repertoire_id"
#' @param count Column with the counts, e.g. "consensus_count",
#'   "duplicate_count" or "umi_count"
#'
#' @returns List with
#'   * sequences: data frame with id, cdr3_aa, v_call and j_call, one row per sequence
#'
#'   * counts: matrix of counts, one row per sequence (rownames "V_CDR3_J") and one column per sample, in the same order as the rows of sequences
#' @noRd
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
  # get levels for factors, order of first appearance, so counts sorting does not depend on locale-specific sorting
  ids <- unique(id)
  samples <- unique(df[[sample]])

  # sum counts per sequence and sample (missing = 0); rowsum() instead of
  # tapply(), which is ~15x slower on large data
  # matrix cell of each input row: column-major position (row + (column - 1) * n_rows)
  cell <- (match(df[[sample]], samples) - 1) * length(ids) + match(id, ids)
  # sum counts per cell (rows of df that are the same sequence in the same sample)
  sums <- rowsum(df[[count]], cell, reorder = FALSE)
  # all other cells (sequence not observed in that sample) stay 0
  counts <- matrix(0, length(ids), length(samples), dimnames = list(ids, samples))
  counts[as.numeric(rownames(sums))] <- sums

  # fix column names
  first <- match(ids, id)
  sequences <- data.frame(id = ids,
                          cdr3_aa = df[[cdr3]][first],
                          v_call = df[[v]][first],
                          j_call = df[[j]][first])

  list(sequences = sequences, counts = counts)
}
