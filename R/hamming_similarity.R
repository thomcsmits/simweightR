hamming_similarity <- function(vj_data) {
  dist.m <- stringdistmatrix(vj_data$cdr3_aa, vj_data$cdr3_aa, method = "hamming")
  rownames(dist.m) <- vj_data$cdr3_aa
  dist.m <- dist.m / nchar(vj_data[1, ]$cdr3_aa)
  sim <- 1 - dist.m
  return(sim)
}
