#' Add column of adjusted counts based on similarity
#'
#' For cdr3 amino acid length in the sample, and for each existing VJ combination,
#' a similarity score matrix is calculated for the respective cdr3s. This similarity
#' matrix is used to generate updated counts via the update_read_counts() function.
#'
#' @param data.1 Dataframe with counts, AIRR format.
#' @inheritParams TCRsimilift_calculate
#'
#' @returns data.1 with added column wrc of adjusted counts.
#' @export
#'
compute_dist_matrix <- function(data.1, sim_method, cutoff=0.8) {

  #Find interval of relevant cdr3 amino acid chain lengths
  min.l <- min(data.1$length)
  max.l <- max(data.1$length)

  data.new <- c()
  # For each cdr3 length in the dataset
  for (seq_length in c(min.l:max.l)) {
    data.l <- data.1[data.1$length == seq_length, ]
    if (nrow(data.l) > 1) {
      unique_vj <- unique(paste(data.l$v_call, data.l$j_call, sep = "_"))
      # For each relevant combination of V and J segments
      for (vj_group in unique_vj) {
        vj_data <- data.l[paste(data.l$v_call, data.l$j_call, sep = "_") == vj_group, ]

        if (nrow(vj_data) > 1) {
          cat("Processing sequences of length", seq_length, "for VJ group", vj_group, "\n")
          # Compute hamming distance-based similarity within the VJ group
          if (sim_method == "HAMMING") {
            sim <- hamming_similarity(vj_data)
          }
          # or BLOSUM-based
          else if (sim_method == "BLOSUM") {
            sim <- blosum_similarity(vj_data)
          }
          else {
            stop("sim_method must be HAMMING or BLOSUM") #throw error if similarity method not valid
          }
          # add new column to dataframe, containing similarity adjusted counts
          vj_data$wrc <- update_read_counts(read_counts = vj_data$consensus_count,
                                            similarity_matrix = sim,
                                            cutoff = cutoff)
        } else {
          vj_data$wrc <- vj_data$consensus_count
        }
        data.new <- base::rbind(data.new, vj_data)
      }
    } else {
      data.l$wrc <- data.l$consensus_count
      data.new <- base::rbind(data.new, data.l)
    }
  }
  return(data.new)
}
