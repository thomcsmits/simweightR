datacheck <- function(df) {
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
}
