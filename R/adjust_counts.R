#' Use simweightR
#'
#' Convenience function to do all simweightR data processing in one call.
#' The function checks the data to ensure format, prepares the data for processing,
#' and runs weight_counts to return an updated dataframe with similarity-adjusted
#' counts.
#'
#' @param df Input dataframe of AIRR formatted immunological data.
#' @param sim_method Similarity method to use. Either \code{"HAMMING"}, \code{"BLOSUM"},
#'   or a custom function. A custom function receives a dataframe of unique sequences
#'   for each VJ+length group (with at least \code{cdr3_aa} and \code{length} columns)
#'   and must return a square numeric matrix of similarity scores between 0 and 1.
#' @param cutoff Minimum similarity score to consider two sequences neighbours. Between 0 and 1. Default 0.8 .
#'
#' @returns Returns dataframe with extra column of adjusted counts based on similarity.
#' @export
#'
#' @examples
#' results <- adjust_counts(mouse_PBSvTCZ_data_minisubset, sim_method="HAMMING", cutoff = 0.77)
#'
adjust_counts <- function(df,
                          sim_method = "HAMMING",
                          cutoff = 0.8) {
  datacheck(df)
  df2 <- dataprep(df)
  weight_counts(df2, sim_method = sim_method, cutoff = cutoff)
}
