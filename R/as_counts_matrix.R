#' Make count matrix from outputs of adjust_counts
#'
#' This function takes a dataframe which is output by adjust_counts,
#' which may contain many columns besides the information about counts,
#' sequence IDs and samples. It then outputs a counts matrix, where each row
#' refers to a sequence, each columns is a sample, and the values therein are
#' counts. Such a count matrix can then be passed on into conventional
#' DGE workflows such as edgeR and deseq2, as well as other methods
#' for analyzing DGE such as the Wilcoxon test.
#'
#' See \link{adjust_counts} for a full example of the DGE workflow.
#'
#' @param new.data A dataframe outputted by \link{adjust_counts}
#' @param doFilter Filter out sequences that do not have at least 2 counts? TRUE/FALSE.
#'
#' @returns Count matrix containing weighted counts for each sequence and sample.
#' @export
#'
#' @examples
#' results <- adjust_counts(mouse_PBSvTCZ_data_minisubset)
#' count_matrix <- as_counts_matrix(results)
#'
#'



as_counts_matrix <- function(new.data, doFilter=FALSE) {
  if (doFilter == TRUE) { #if filtering, we need to process unweighted data too
    ## Save unweighted counts
    new.rc.unwghtd <- new.data[, c("sequence_id", "consensus_count", "sample_processing_id")]
    new.rc.unwghtd <- data.frame(new.rc.unwghtd)
    wide_df <- reshape(new.rc.unwghtd,
                       idvar = "sequence_id",
                       timevar = "sample_processing_id",
                       direction = "wide")
    rownames(wide_df) <- wide_df$sequence_id
    wide_df <- wide_df[, -1]
    # Set the column names
    colnames(wide_df) <- unique(new.rc.unwghtd$sample_processing_id)
    # Convert to matrix
    result_matrix_unweighted <- as.matrix(wide_df)
    result_matrix_unweighted[is.na(result_matrix_unweighted)] <- 0
    # Filter rows that occur in at least 2 columns
    valid_rows <- rowSums(result_matrix_unweighted > 0) >= 2#ceiling(sample_size/4)
  }

  ## Save weighted counts
  new.rc <- new.data[, c("sequence_id", "wrc", "sample_processing_id")]
  new.rc <- data.frame(new.rc)
  wide_df <- reshape(new.rc,
                     idvar = "sequence_id",
                     timevar = "sample_processing_id",
                     direction = "wide")
  rownames(wide_df) <- wide_df$sequence_id
  wide_df <- wide_df[, -1]
  # Set the column names
  colnames(wide_df) <- unique(new.rc$sample_processing_id)
  # Convert to matrix
  result_matrix_weighted <- as.matrix(wide_df)
  result_matrix_weighted[is.na(result_matrix_weighted)] <- 0
  # Use the same valid rows to filter the weighted matrix
  if (doFilter == TRUE) {
    return(result_matrix_weighted[valid_rows, ])
  }
  else {
    return(result_matrix_weighted)
  }
}
