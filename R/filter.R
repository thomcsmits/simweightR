#' Keep sequences observed in enough samples
#'
#' Keeps sequences with a count of at least min_count in at least min_samples
#' samples, judged on the unweighted counts. counts, sequences and (if present)
#' weighted are subset together, so they stay aligned. Can be used before or
#' after weight_counts().
#'
#' @param x List with counts and sequences (and optionally weighted), e.g. from
#'   prepare_counts() or weight_counts()
#' @param min_samples Minimum number of samples a sequence must be observed in
#' @param min_count Minimum count for a sequence to be observed in a sample
#'
#' @returns x with only the kept sequences
#'
#' @export
filter_counts <- function(x, min_samples = 2, min_count = 1) {
  # judged on the unweighted counts
  keep <- rowSums(x$counts >= min_count) >= min_samples

  x$counts <- x$counts[keep, , drop = FALSE]
  x$sequences <- x$sequences[keep, , drop = FALSE]
  if (!is.null(x$weighted)) {
    x$weighted <- x$weighted[keep, , drop = FALSE]
  }
  x
}
