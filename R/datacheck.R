#' Check if dataset has all necessary columns
#'
#' This internal function of the package verifies that the input dataframe contains the
#' needed AIRR columns. It does not verify their contents, only their existence.
#'
#' The function will only have outputs if there are any errors or warnings, otherwise
#' the dataset is considered OK.
#'
#' See \link{adjust_counts} for a full example of the DGE workflow.
#'
#' @param df Dataframe of immunological data.
#'
#' @export
#'
datacheck <- function(df) {

  print("Checking dataset...")

  necessary_cols <- c("cdr3_aa",
                      "consensus_count",
                      "v_call",
                      "j_call",
                      "sample_processing_id",
                      "sequence_id")

  if (!("cdr3" %in% colnames(df))) {
    print("Note: cdr3 nucleotide sequence not present in dataset. Data processing can
          proceed correctly regardless.")
  }

  if (all(necessary_cols %in% colnames(df)) == FALSE) {
    stop("Dataset incomplete. The dataset must have at least the following columns:
         cdr3_aa, consensus_count, v_call, j_call, sample_processing_id, sequence_id")
  }
  else {
    print("Dataset contains all necessary columns.")
  }
}
