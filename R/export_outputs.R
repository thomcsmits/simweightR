#' Export results of TCRsimilift
#'
#' Creates output directory if needed, and saves in it R datasets corresponding
#' to the original and respectively to the adjusted counts.
#'
#' @param new.data Dataframe produced by the TCRsimilift_calculate function.
#' @param output_directory Name of output folder. "outputs" by default.
#'
#' @export
#'
export_outputs <- function(new.data, output_directory = "outputs") {

  #If output directory does not exist, create it
  if (!(dir.exists(output_directory))) {
    dir.create(output_directory)
  }

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

  write_rds(result_matrix, file=paste0(output_directory,
                                       "/matrix-unweighted.Rds"))

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

  write_rds(result_matrix, file=paste0(output_directory,
                                       "/matrix-weighted-average.Rds"))

  ## Save csv file with full outputs of data processing
  write.csv(new.data,
            gzfile(paste0(output_directory, "/update_data_TCRsimilift.csv.gz")))
}
