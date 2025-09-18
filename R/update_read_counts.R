update_read_counts <- function(read_counts, similarity_matrix) {
  updated_read_counts <- read_counts
  for (i in 1:length(read_counts)) {
    similar_indices <- which(similarity_matrix[i, ] >= 0.8 ) #HARD-CODED VALUE!
    weighted_sum <- 0
    total_weight <- 0

    if (read_counts[i] > 0) {
      weighted_sum <- read_counts[i] + sum(similarity_matrix[i, similar_indices] * (read_counts[similar_indices]))
      total_weight <- 1 + sum(similarity_matrix[i, similar_indices])
    } else {
      weighted_sum <- sum(similarity_matrix[i, similar_indices] * read_counts[similar_indices])
      total_weight <- sum(similarity_matrix[i, similar_indices])
    }
    if (total_weight > 0) {
      updated_read_counts[i] <- weighted_sum / total_weight
    }
  }
  return(updated_read_counts)
}
