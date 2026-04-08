#' Add column of adjusted counts based on similarity
#'
#' For each cdr3 amino acid length and V-J combination, a similarity matrix is
#' calculated once for all unique sequences in that group, then applied to each
#' sample. This avoids recomputing the similarity matrix per sample.
#'
#' See \link{TCRsimilift_calculate} for a full example of the DGE workflow.
#'
#' @param data.1 Dataframe with counts for all samples, AIRR format.
#' @inheritParams TCRsimilift_calculate
#'
#' @returns data.1 with added column wrc of adjusted counts.
#'
compute_dist_matrix <- function(data.1, sim_method, cutoff=0.8) {

  min.l <- min(data.1$length)
  max.l <- max(data.1$length)

  data.new <- c()
  for (seq_length in seq(min.l, max.l)) {
    data.l <- data.1[data.1$length == seq_length, ]
    unique_vj <- unique(paste(data.l$v_call, data.l$j_call, sep = "_"))

    for (vj_group in unique_vj) {
      vj_data <- data.l[paste(data.l$v_call, data.l$j_call, sep = "_") == vj_group, ]
      unique_seqs <- unique(vj_data$cdr3_aa)

      if (length(unique_seqs) > 1) {
        cat("Processing sequences of length", seq_length, "for VJ group", vj_group, "\n")

        # Compute similarity matrix ONCE for all unique sequences in this group
        seq_df <- vj_data[match(unique_seqs, vj_data$cdr3_aa), ]
        if (sim_method == "HAMMING") {
          sim <- hamming_similarity(seq_df)
        } else if (sim_method == "BLOSUM") {
          sim <- blosum_similarity(seq_df)
        } else {
          stop("sim_method must be HAMMING or BLOSUM")
        }

        # Apply the same similarity matrix to each sample's counts
        for (s in unique(vj_data$sample_processing_id)) {
          sample_data <- vj_data[vj_data$sample_processing_id == s, ]
          ord <- match(sample_data$cdr3_aa, unique_seqs)
          sample_data$wrc <- update_read_counts(sample_data$consensus_count,
                                                sim[ord, ord, drop = FALSE],
                                                cutoff)
          data.new <- base::rbind(data.new, sample_data)
        }
      } else {
        vj_data$wrc <- vj_data$consensus_count
        data.new <- base::rbind(data.new, vj_data)
      }
    }
  }
  return(data.new)
}
