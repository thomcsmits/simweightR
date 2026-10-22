#' Similarity based weighting of counts data
#'
#' simweight prepares and weights counts in one call
#'
#' Runs [prepare_counts()] and [weight_counts()]
#'
#' @param df AIRR formatted data frame with one row per sequence
#' @param threshold Minimum similarity for two sequences to be neighbours,
#'   between 0 and 1
#' @param method Similarity method: "hamming" (default), "blosum", or a
#'   function that takes a character vector of CDR3 sequences of one group
#'   (all the same length) and returns their n x n similarity matrix, with
#'   values between 0 and 1, in the same order
#' @param nthread Number of threads used by stringdist; only used by method
#'   "hamming"
#' @param n_cores Number of V/J/length groups processed in parallel, with
#'   parallel::mclapply(). Not available on Windows: there it falls back to 1
#'   with a warning. The total number of threads is n_cores x nthread; on a busy
#'   shared machine more threads can be slower. Warnings from the workers are
#'   not shown.
#' @param cdr3 Column with the CDR3 amino acid sequence, e.g. "cdr3_aa" or
#'   "junction_aa" (junction = CDR3 plus the conserved flanking C and F/W)
#' @param v Column with the V gene
#' @param j Column with the J gene
#' @param sample Column with the sample identifier, e.g.
#'   "sample_processing_id" or "repertoire_id"
#' @param count Column with the counts, e.g. "consensus_count",
#'   "duplicate_count" or "umi_count"
#'
#' @returns List with
#'   * sequences: data frame with id, cdr3_aa, v_call and j_call, one row per
#'     sequence
#'
#'   * counts: matrix of counts, one row per sequence (rownames "V_CDR3_J") and
#'     one column per sample, in the same order as the rows of sequences
#'
#'   * weighted: matrix of similarity-weighted counts, same dimensions and row
#'     order as counts
#'
#' @seealso [prepare_counts()] and [weight_counts()] to run the two steps
#'   separately, e.g. to try several thresholds without repeating the
#'   preparation.
#'
#' @export
simweight <- function(df,
                      threshold = 0.8,
                      method = "hamming",
                      nthread = 1,
                      n_cores = 1,
                      cdr3 = "cdr3_aa",
                      v = "v_call",
                      j = "j_call",
                      sample = "sample_processing_id",
                      count = "consensus_count") {
  x <- prepare_counts(df, cdr3 = cdr3, v = v, j = j, sample = sample,
                      count = count)
  x <- weight_counts(x, threshold = threshold, method = method,
                     nthread = nthread, n_cores = n_cores)
  x
}
