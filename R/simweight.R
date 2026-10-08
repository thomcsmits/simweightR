#' Prepare and weight counts in one call
#'
#' Runs prepare_counts() and weight_counts(). Filtering is not included yet
#' (whether to filter before or after weighting is still open); call
#' filter_counts() on the result.
#'
#' @param df AIRR formatted data frame, see prepare_counts()
#' @param threshold Minimum similarity for two sequences to be neighbours,
#'   between 0 and 1
#' @param method Similarity method: "hamming" (default), "blosum", or a
#'   function, see weight_counts()
#' @param nthread Number of threads for stringdist, see hamming_similarity();
#'   only used by method "hamming"
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
#' @returns List with counts, sequences and weighted, see weight_counts()
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
  x <- prepare_counts(df, cdr3 = cdr3, v = v, j = j, sample = sample, count = count)
  x <- weight_counts(x, threshold = threshold, method = method, nthread = nthread, n_cores = n_cores)
  # filter tbd
  # x <- filter_counts(x, min_samples = min_samples, min_count = min_count)
  x
}