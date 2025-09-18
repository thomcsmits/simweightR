compute_dist_matrix <- function(data.1, z) {
  cat("Processing mouse", z, "\n")

  # Parse V gene and J gene from seqid
  # data.1$v_gene <- gsub("_.*", "", data.1$seqid)    #get vgene jgene directly instead!
  # data.1$j_gene <- gsub(".*_", "", data.1$seqid)

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
          # Compute hamming distance within the VJ group
          dist.m <- stringdistmatrix(vj_data$cdr3_aa, vj_data$cdr3_aa, method = "hamming")
          rownames(dist.m) <- vj_data$cdr3_aa
          dist.m <- dist.m / nchar(vj_data[1, ]$cdr3_aa)
          sim <- 1 - dist.m
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
