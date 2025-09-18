export_outputs <- function(new.data) {
  ## Save unweighted counts
  new.rc <- new.data[,c("sequence_id", "consensus_count", "sample_processing_id")]
  new.rc <- data.frame(new.rc)
  colnames(new.rc) <- c("junction_aa", "counts", "i")

  wide_df <- reshape(new.rc, idvar = "junction_aa", timevar = "i", direction = "wide")
  rownames(wide_df) <- wide_df$junction_aa

  wide_df <- wide_df[, -1]

  # Set the column names
  colnames(wide_df) <- unique(new.rc$i)

  # Convert to matrix
  result_matrix <- as.matrix(wide_df)
  result_matrix[is.na(result_matrix)] <- 0

  sample_size = 2

  write_rds(result_matrix, file=paste0("data/matrix-unweighted-", sample_size, ".Rds"))

  ## Save weighted counts
  new.rc <- new.data[,c("sequence_id", "wrc", "sample_processing_id")]
  new.rc <- data.frame(new.rc)
  colnames(new.rc) <- c("junction_aa", "counts", "i")

  wide_df <- reshape(new.rc, idvar = "junction_aa", timevar = "i", direction = "wide")
  rownames(wide_df) <- wide_df$junction_aa

  wide_df <- wide_df[, -1]

  # Set the column names
  colnames(wide_df) <- unique(new.rc$i)

  # Convert to matrix
  result_matrix <- as.matrix(wide_df)
  result_matrix[is.na(result_matrix)] <- 0

  write_rds(result_matrix, file=paste0("data/matrix-weighted-average-", sample_size, ".Rds"))

  return(":)")
}
