blosum_similarity <- function(vj_data) {
  seq1 <- expand.grid(vj_data$cdr3_aa, vj_data$cdr3_aa)
  dist.m <- matrix(score(pairwiseAlignment(seq1$Var1, seq1$Var2,
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
