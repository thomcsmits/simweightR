#' Calculate similarity matrix for input sequences based on BLOSUM distance
#'
#' Calculates BLOSUM62-based distance matrix for the given sequences, which is then
#' employed to generate a similarity matrix, which is returned.
#'
#' See \link{TCRsimilift_calculate} for a full example of the DGE workflow.
#'
#' @param vj_data Data subset containing all sequences with same length, same V and J calls.
#'
#' @returns Similarity score matrix.
#'
blosum_similarity <- function(vj_data) {

  #Check whether data is appropriate
  if (!(nrow(vj_data) > 1)) {
    stop("At least two sequences needed to calculate similarity.")
  }
  if (!(length(unique(vj_data$length)) == 1)) {
    stop("All amino acid sequences must have same length to calculate similarity.")
  }

  seq1 <- expand.grid(vj_data$cdr3_aa, vj_data$cdr3_aa)
  dist.m <- matrix(BiocGenerics::score(Biostrings::pairwiseAlignment(
                                           seq1$Var1, seq1$Var2,
                                           substitutionMatrix = "BLOSUM62",
                                           gapOpening = -2,
                                           gapExtension = -8,
                                           scoreOnly = FALSE)),
                   nrow = length(vj_data$cdr3_aa))
  rownames(dist.m) <- vj_data$cdr3_aa
  colnames(dist.m) <- vj_data$cdr3_aa
  dist.m <- dist.m / nchar(vj_data[1, ]$cdr3_aa)
  sim <- apply(dist.m, 2, function(col) (col - min(col)) / (max(col) - min(col)))
  return(sim)
}
