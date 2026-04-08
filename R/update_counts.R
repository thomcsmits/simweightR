#' Update read counts based on similarity
#'
#' Given the input read counts and similarity matrix for them, updated counts are
#' generated based on the weighted-sum rule.
#'
#' See \link{adjust_counts} for a full example of the DGE workflow.
#'
#' @param read_counts A dataframe with read counts for each sequence ID
#' @param similarity_matrix A matrix of similarity scores for each sequence ID pair
#' @inheritParams adjust_counts
#'
#' @returns A dataframe.
#'
update_counts <- function(read_counts, similarity_matrix, cutoff=0.8) {
  masked_sim <- Matrix::Matrix(similarity_matrix * (similarity_matrix >= cutoff), sparse = TRUE)
  neighbor_sum <- as.vector(masked_sim %*% read_counts)
  neighbor_weight <- Matrix::rowSums(masked_sim)

  has_counts <- read_counts > 0
  numerator <- ifelse(has_counts, read_counts + neighbor_sum, neighbor_sum)
  denominator <- ifelse(has_counts, 1 + neighbor_weight, neighbor_weight)

  ifelse(denominator > 0, numerator / denominator, read_counts)
}
