compute_dist_matrix <- function(data.1, z, sim_method) {

  min.l <- min(data.1$length)
  max.l <- max(data.1$length)

  data.new <- c()
  for (seq_length in c(min.l:max.l)) {
    data.l <- data.1[data.1$length == seq_length, ]
    if (nrow(data.l) > 1) {
      unique_vj <- unique(paste(data.l$v_call, data.l$j_call, sep = "_"))
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
          vj_data$wrc <- update_read_counts(vj_data$consensus_count, sim)
        } else {
          vj_data$wrc <- vj_data$consensus_count
        }
        data.new <- rbind(data.new, vj_data)
      }
    } else {
      data.l$wrc <- data.l$consensus_count
      data.new <- rbind(data.new, data.l)
    }
  }
  return(data.new)
}
