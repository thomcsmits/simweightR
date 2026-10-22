#' Similarity-weighted counts for one V/J/length group
#'
#' For sequence s_i, the weighted count is C'_i = sum_j M'_ij C_j / sum_j M'_ij,
#' where M'_ij = M_ij if M_ij >= threshold, else 0. Self-similarity M'_ii = 1.
#'
#' @param counts Matrix of counts, one row per sequence and one column per
#'   sample
#' @param similarity Square matrix of pairwise similarities, rows in the same
#'   order as counts, e.g. from hamming_similarity()
#' @param threshold Minimum similarity for two sequences to be neighbours,
#'   between 0 and 1
#'
#' @returns Matrix of weighted counts, same dimensions as counts
#'
#' @noRd
weight_group <- function(counts, similarity, threshold) {
  if (!is.matrix(counts)) {
    stop("counts must be a matrix with one row per sequence and one column per sample.")
  }
  if (nrow(similarity) != nrow(counts)) {
    stop("similarity must have one row per sequence in counts.")
  }
  if (threshold < 0 || threshold > 1) {
    stop("threshold must be between 0 and 1.")
  }

  # neighbours: similarity >= threshold, with a small tolerance so that pairs
  # 1e-10 is for floating point errors
  # (e.g. 1 - 3/20 vs seq(0.5, 0.95, by = 0.05)[8])
  neighbours <- similarity * (similarity >= threshold - 1e-10)

  # the sequence itself is counted once, with weight 1
  if (any(abs(diag(similarity) - 1) > 1e-10)) {
    warning("similarity does not have 1 on the diagonal (self-similarity); setting it to 1")
    diag(neighbours) <- 1
  }

  (neighbours %*% counts) / rowSums(neighbours)
}

#' Similarity-weighted counts for all sequences
#'
#' Splits the sequences into groups with the same V gene, J gene and CDR3
#' length, computes the similarity within each group and replaces the counts by
#' similarity-weighted counts.
#'
#' @param x List with counts (matrix, one row per sequence and one column per
#'   sample) and sequences (data frame with id, cdr3_aa, v_call and j_call, in
#'   the same order as the rows of counts), e.g. from prepare_counts()
#' @inheritParams simweight
#'
#' @returns x with an extra element weighted: matrix of weighted counts, same
#'   dimensions and row order as x$counts
#'
#' @examples
#' x <- prepare_counts(mouse_pbs_srbc_mini)
#'
#' # try several thresholds without repeating the preparation
#' x1 <- weight_counts(x, threshold = 0.8)
#' x2 <- weight_counts(x, threshold = 0.9)
#' head(x2$weighted[, 1:4])
#'
#' @export
weight_counts <- function(x, threshold, method = "hamming", nthread = 1,
                          n_cores = 1) {
  if (!identical(rownames(x$counts), x$sequences$id)) {
    stop("rows of x$counts and x$sequences are not in the same order.")
  }

  # similarity function per group: built in by name, or given by the user;
  # each is called with only the sequences of one group
  if (is.function(method)) {
    similarity <- method
  } else if (identical(method, "hamming")) {
    similarity <- function(sequences) {
      hamming_similarity(sequences, nthread = nthread)
    }
  } else if (identical(method, "blosum")) {
    similarity <- function(sequences) blosum_similarity(sequences)
  } else {
    stop('method must be "hamming", "blosum" or a function.')
  }

  # row numbers per V/J/length group; drop = TRUE skips empty combinations
  groups <- split(seq_len(nrow(x$counts)),
                  list(x$sequences$v_call, x$sequences$j_call,
                       nchar(x$sequences$cdr3_aa)),
                  drop = TRUE)

  # weighted counts of one group
  weight_one <- function(rows) {
    sim <- similarity(x$sequences$cdr3_aa[rows])
    check_similarity(sim, length(rows))
    weight_group(x$counts[rows, , drop = FALSE], sim, threshold)
  }

  # mclapply() cannot fork on Windows
  if (n_cores > 1 && .Platform$OS.type == "windows") {
    warning("n_cores > 1 is not supported on Windows; using n_cores = 1.")
    n_cores <- 1
  }

  if (n_cores > 1) {
    results <- parallel::mclapply(groups, weight_one, mc.cores = n_cores,
                                  mc.preschedule = FALSE)
    # errors in a worker come back as try-error objects
    failed <- vapply(results, inherits, logical(1), what = "try-error")
    if (any(failed)) {
      stop(results[[which(failed)[1]]])
    }
  } else {
    results <- lapply(groups, weight_one)
  }

  weighted <- x$counts
  for (k in seq_along(groups)) {
    weighted[groups[[k]], ] <- results[[k]]
  }

  # weighted must have the same rows, in the same order, as counts
  if (!identical(dim(x$counts), dim(weighted)) ||
        !identical(rownames(x$counts), rownames(weighted))) {
    stop("weighted counts do not match the rows of x$counts.")
  }

  x$weighted <- weighted
  x
}
