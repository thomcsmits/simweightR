#' Calculate similarity matrix for input sequences based on Hamming distance
#'
#' See \link{TCRsimilift_calculate} for a full example of the DGE workflow.
#'
#' @param vj_data Data subset containing all sequences with same length, same V and J calls.
#'
#' @returns Similarity score matrix.
#'
hamming_similarity <- function(vj_data) {

  #Check whether data is appropriate
  if (!(nrow(vj_data) > 1)) {
    stop("At least two sequences needed to calculate similarity.")
  }
  if (!(length(unique(vj_data$length)) == 1)) {
    stop("All amino acid sequences must have same length to calculate similarity.")
  }

  dist.m <- stringdist::stringdistmatrix(vj_data$cdr3_aa,
                                         vj_data$cdr3_aa,
                                         method = "hamming")
  rownames(dist.m) <- vj_data$cdr3_aa
  dist.m <- dist.m / nchar(vj_data[1, ]$cdr3_aa)
  sim <- 1 - dist.m
  return(sim)
}
