# Regression tests: verify the package produces the same wrc values as the
# original analysis code in dge-analysis/figures/figure-6-new/code/util.R.
#
# Column name mapping (package -> figure-6):
#   cdr3_aa              -> junction_aa
#   sequence_id          -> seqid  (already encoded as v_call_cdr3_aa_j_call)
#   consensus_count      -> counts
#   sample_processing_id -> i

# ---------------------------------------------------------------------------
# Inline re-implementation of figure-6 util.R (base R, no tidyverse)
# ---------------------------------------------------------------------------

fig6_update_read_counts <- function(read_counts, similarity_matrix, threshold) {
  masked_sim      <- similarity_matrix * (similarity_matrix >= threshold)
  neighbor_sum    <- as.vector(masked_sim %*% read_counts)
  neighbor_weight <- rowSums(masked_sim)

  has_counts  <- read_counts > 0
  numerator   <- ifelse(has_counts, read_counts + neighbor_sum, neighbor_sum)
  denominator <- ifelse(has_counts, 1 + neighbor_weight, neighbor_weight)

  ifelse(denominator > 0, numerator / denominator, read_counts)
}

fig6_sim_func_hamming <- function(junction_aa, seq_length) {
  n   <- length(junction_aa)
  idx <- which(upper.tri(matrix(0, n, n)), arr.ind = TRUE)
  d   <- stringdist::stringdist(junction_aa[idx[, 1]],
                                junction_aa[idx[, 2]],
                                method = "hamming") / seq_length
  dist.m <- matrix(0, n, n)
  dist.m[idx]                          <- d
  dist.m[idx[, c(2, 1), drop = FALSE]] <- d
  rownames(dist.m) <- junction_aa
  1 - dist.m
}

# Run a figure-6-style pipeline on package-format data.
fig6_pipeline <- function(df, threshold = 0.8) {
  df$length <- nchar(df$cdr3_aa)
  data <- stats::aggregate(consensus_count ~ sequence_id + sample_processing_id,
                           df, FUN = sum, drop = FALSE)
  data$consensus_count[is.na(data$consensus_count)] <- 0

  meta <- unique(df[c("sequence_id", "length", "v_call", "j_call", "cdr3_aa")])
  data$length  <- meta$length [match(data$sequence_id, meta$sequence_id)]
  data$v_call  <- meta$v_call [match(data$sequence_id, meta$sequence_id)]
  data$j_call  <- meta$j_call [match(data$sequence_id, meta$sequence_id)]
  data$cdr3_aa <- meta$cdr3_aa[match(data$sequence_id, meta$sequence_id)]
  data$vj_group <- paste(data$v_call, data$j_call, sep = "_")

  samples  <- unique(data$sample_processing_id)
  data$wrc <- NA_real_

  for (seq_length in unique(data$length)) {
    for (vj in unique(data$vj_group[data$length == seq_length])) {
      rows_vj     <- which(data$length == seq_length & data$vj_group == vj)
      unique_seqs <- unique(data$cdr3_aa[rows_vj])

      if (length(unique_seqs) <= 1) {
        data$wrc[rows_vj] <- data$consensus_count[rows_vj]
        next
      }

      sim <- fig6_sim_func_hamming(unique_seqs, seq_length)

      for (s in samples) {
        rows_s <- rows_vj[data$sample_processing_id[rows_vj] == s]
        ord    <- match(data$cdr3_aa[rows_s], unique_seqs)
        data$wrc[rows_s] <- fig6_update_read_counts(
          data$consensus_count[rows_s], sim[ord, ord, drop = FALSE], threshold)
      }
    }
  }

  data$vj_group <- NULL
  data
}

# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

test_that("fig6_sim_func_hamming diagonal is 1 for n=2 sequences", {
  seqs <- c("ASGDKRGEQY", "ASGDSRNEQY")
  m    <- fig6_sim_func_hamming(seqs, 10)
  expect_equal(diag(m), c(1, 1))
})

test_that("hamming_similarity (package) and fig6_sim_func_hamming give identical matrices", {
  vj <- data.frame(
    cdr3_aa = c("ASGDKRGEQY", "ASGDSRNEQY", "ASGDTRGEQY"),
    length = 10,
    v_call = "TRBV13-2",
    j_call = "TRBJ2-7"
  )
  pkg_sim <- hamming_similarity(vj)
  fig_sim <- fig6_sim_func_hamming(vj$cdr3_aa, 10)

  expect_equal(unname(pkg_sim), unname(fig_sim))
})

test_that("update_read_counts and fig6_update_read_counts produce identical results", {
  counts <- c(10, 5, 0, 3)
  sim    <- matrix(c(
    1.0, 0.9, 0.5, 0.3,
    0.9, 1.0, 0.6, 0.4,
    0.5, 0.6, 1.0, 0.8,
    0.3, 0.4, 0.8, 1.0
  ), nrow = 4)

  pkg_result <- update_counts(counts, sim, cutoff = 0.8)
  fig_result <- fig6_update_read_counts(counts, sim, threshold = 0.8)

  expect_equal(pkg_result, fig_result)
})

test_that("full pipeline: package and figure-6 logic return identical wrc values", {
  pkg_result <- adjust_counts(mouse_PBSvTCZ_data_minisubset,
                                      sim_method = "HAMMING",
                                      cutoff = 0.8)
  fig_result <- fig6_pipeline(mouse_PBSvTCZ_data_minisubset, threshold = 0.8)

  key_pkg <- paste(pkg_result$sequence_id, pkg_result$sample_processing_id)
  key_fig <- paste(fig_result$sequence_id, fig_result$sample_processing_id)
  fig_wrc <- fig_result$wrc[match(key_pkg, key_fig)]

  expect_equal(pkg_result$wrc, fig_wrc, tolerance = 1e-9)
})
