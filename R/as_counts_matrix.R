#' Make count matrix from outputs of adjust_counts
#'
#' This function takes a dataframe which is output by adjust_counts,
#' which may contain many columns besides the information about counts,
#' sequence IDs and samples. It then outputs a counts matrix, where each row
#' refers to a sequence, each column is a sample, and the values therein are
#' counts. Such a count matrix can then be passed on into conventional
#' DGE workflows such as edgeR and deseq2, as well as other methods
#' for analyzing DGE such as the Wilcoxon test.
#'
#' See \link{adjust_counts} for a full example of the DGE workflow.
#'
#' @param new.data A dataframe outputted by \link{adjust_counts}
#' @param doFilter Filter out sequences that do not have at least 2 counts? TRUE/FALSE.
#' @param weighted Return similarity-adjusted counts (\code{wrc})? If \code{FALSE},
#'   returns original counts (\code{consensus_count}). Default \code{TRUE}.
#'
#' @returns Count matrix containing counts for each sequence and sample.
#' @export
#'
#' @examples
#' results <- adjust_counts(mouse_PBSvTCZ_data_minisubset)
#' count_matrix <- as_counts_matrix(results)
#' count_matrix_unweighted <- as_counts_matrix(results, weighted = FALSE)
#'


as_counts_matrix <- function(new.data, doFilter = FALSE, weighted = TRUE) {
  count_col <- if (weighted) "wrc" else "consensus_count"

  if (doFilter) {
    ## Build unweighted matrix for filtering
    new.rc.unwghtd <- new.data[, c("sequence_id", "consensus_count", "sample_processing_id")]
    new.rc.unwghtd <- data.frame(new.rc.unwghtd)
    wide_df <- reshape(new.rc.unwghtd,
                       idvar = "sequence_id",
                       timevar = "sample_processing_id",
                       direction = "wide")
    rownames(wide_df) <- wide_df$sequence_id
    wide_df <- wide_df[, -1]
    colnames(wide_df) <- unique(new.rc.unwghtd$sample_processing_id)
    result_matrix_unweighted <- as.matrix(wide_df)
    result_matrix_unweighted[is.na(result_matrix_unweighted)] <- 0
    valid_rows <- rowSums(result_matrix_unweighted > 0) >= 2
  }

  ## Build requested count matrix
  new.rc <- new.data[, c("sequence_id", count_col, "sample_processing_id")]
  new.rc <- data.frame(new.rc)
  wide_df <- reshape(new.rc,
                     idvar = "sequence_id",
                     timevar = "sample_processing_id",
                     direction = "wide")
  rownames(wide_df) <- wide_df$sequence_id
  wide_df <- wide_df[, -1]
  colnames(wide_df) <- unique(new.rc$sample_processing_id)
  result_matrix <- as.matrix(wide_df)
  result_matrix[is.na(result_matrix)] <- 0

  if (doFilter) {
    return(result_matrix[valid_rows, ])
  }
  result_matrix
}
