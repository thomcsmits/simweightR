
#' Hamming similarity between sequences of equal length
#'
#' For sequences s_i and s_j of length n, the similarity is 1 - h_ij / n,
#' where h_ij is the Hamming distance.
#'
#' @param sequences Vector with sequences of same length
#' @param nthread Number of threads for stringdist. Default 1; the stringdist
#'   default (most cores) can be much slower on a busy shared machine.
#'
#' @returns Hamming similarity in matrix of pairwise similarities (diagonal=1)
#'
#' @noRd
hamming_similarity <- function(sequences, nthread = 1) {
  # check whether data is appropriate
  if (!is.character(sequences)) {
    stop("sequences must be a character vector, e.g. vj_data$cdr3_aa.")
  }

  if (!(length(unique(nchar(sequences))) == 1)) {
    stop("All amino acid sequences must have same length to calculate similarity.")
  }

  1 - stringdist::stringdistmatrix(sequences, sequences, method = "hamming",
                                   nthread = nthread) / nchar(sequences[1])
}

#' BLOSUM similarity between sequences of equal length
#'
#' For sequences s_i and s_j of length n, the similarity is
#' `1 - sum_p d(s_i[p], s_j[p]) / n`, where d is the distance between the amino
#' acids at position p (blosum62_distance). Same form as the Hamming
#' similarity, where d is 0 for equal and 1 for different amino acids.
#'
#' @param sequences Vector with sequences of same length
#' @param distance Symmetric amino acid distance matrix with amino acids as row
#'   and column names. Default blosum62_distance.
#'
#' @returns BLOSUM similarity in matrix of pairwise similarities (diagonal=1)
#'
#' @noRd
blosum_similarity <- function(sequences,
                              distance = simweightR::blosum62_distance) {
  # check whether data is appropriate
  if (!is.character(sequences)) {
    stop("sequences must be a character vector, e.g. vj_data$cdr3_aa.")
  }

  if (!(length(unique(nchar(sequences))) == 1)) {
    stop("All amino acid sequences must have same length to calculate similarity.")
  }

  # amino acid at each position as row number in distance: one row per
  # sequence, one column per position
  n_char <- nchar(sequences[1])
  aa <- matrix(match(unlist(strsplit(sequences, "")), rownames(distance)),
               nrow = length(sequences), byrow = TRUE)
  if (anyNA(aa)) {
    stop("sequences contain characters that are not in the distance matrix.")
  }

  # sum the per-position distances over all positions
  d <- unname(distance)
  total <- matrix(0, length(sequences), length(sequences))
  for (p in seq_len(n_char)) {
    total <- total + d[aa[, p], aa[, p], drop = FALSE]
  }

  1 - total / n_char
}

#' Check a similarity matrix returned by a similarity method
#'
#' Every similarity method (hamming_similarity(), blosum_similarity() or a
#' custom function) must return an n x n matrix with values between 0 and 1.
#'
#' @param sim Similarity matrix
#' @param n Number of sequences it should be for
#'
#' @noRd
check_similarity <- function(sim, n) {
  if (!is.matrix(sim)) {
    stop("similarity method must return a matrix.")
  }
  if (nrow(sim) != n || ncol(sim) != n) {
    stop(sprintf("similarity matrix must be %d x %d, got %d x %d.",
                 n, n, nrow(sim), ncol(sim)))
  }
  if (any(sim < 0 | sim > 1)) {
    stop("similarity matrix values must be between 0 and 1.")
  }
}
