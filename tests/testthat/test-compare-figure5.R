# Regression tests: verify the package matches the figure-5 analysis code in
# dge-analysis/figures/figure-5-new/code-dge/prepare-dge-hamming.R.
#
# Figure-5 uses a loop-based per-sample approach that is structurally identical
# to the package. Column name mapping:
#   cdr3_aa              -> junction_aa
#   consensus_count      -> Read.count
#   sequence_id          -> seqid  (v_call_cdr3_aa_j_call)
#   sample_processing_id -> i
#   v_call               -> v_gene  (parsed from seqid via gsub)
#   j_call               -> j_gene  (parsed from seqid via gsub)
#
# No known bugs. The algorithms are line-for-line equivalent, so all wrc
# values should match to floating-point precision.

# ---------------------------------------------------------------------------
# Inline re-implementation of figure-5 prepare-dge-hamming.R (base R only)
# ---------------------------------------------------------------------------

fig5_update_read_counts <- function(read_counts, similarity_matrix, threshold) {
  updated_read_counts <- read_counts
  for (i in seq_along(read_counts)) {
    similar_indices <- which(similarity_matrix[i, ] >= threshold)
    weighted_sum <- 0
    total_weight <- 0

    if (read_counts[i] > 0) {
      weighted_sum <- read_counts[i] +
        sum(similarity_matrix[i, similar_indices] * read_counts[similar_indices])
      total_weight <- 1 + sum(similarity_matrix[i, similar_indices])
    } else {
      weighted_sum <- sum(similarity_matrix[i, similar_indices] * read_counts[similar_indices])
      total_weight <- sum(similarity_matrix[i, similar_indices])
    }

    if (total_weight > 0) {
      updated_read_counts[i] <- weighted_sum / total_weight
    }
  }
  updated_read_counts
}

fig5_compute_hamming_by_length <- function(data.1, threshold = 0.8) {
  data.1$v_gene <- gsub("_.*", "", data.1$seqid)
  data.1$j_gene <- gsub(".*_", "", data.1$seqid)

  min.l <- min(data.1$length)
  max.l <- max(data.1$length)

  data.new <- NULL
  for (seq_length in seq(min.l, max.l)) {
    data.l <- data.1[data.1$length == seq_length, ]
    if (nrow(data.l) > 1) {
      unique_vj <- unique(paste(data.l$v_gene, data.l$j_gene, sep = "_"))
      for (vj_group in unique_vj) {
        vj_data <- data.l[paste(data.l$v_gene, data.l$j_gene, sep = "_") == vj_group, ]
        if (nrow(vj_data) > 1) {
          dist.m <- stringdist::stringdistmatrix(vj_data$junction_aa,
                                                 vj_data$junction_aa,
                                                 method = "hamming")
          rownames(dist.m) <- vj_data$junction_aa
          dist.m <- dist.m / nchar(vj_data[1, ]$junction_aa)
          sim <- 1 - dist.m
          vj_data$wrc <- fig5_update_read_counts(vj_data$Read.count, sim, threshold)
        } else {
          vj_data$wrc <- vj_data$Read.count
        }
        data.new <- rbind(data.new, vj_data)
      }
    } else {
      data.l$wrc <- data.l$Read.count
      data.new <- rbind(data.new, data.l)
    }
  }
  data.new
}

# Run the figure-5 pipeline on package-format data.
fig5_pipeline <- function(df, threshold = 0.8) {
  df$length <- nchar(df$cdr3_aa)
  data <- stats::aggregate(consensus_count ~ sequence_id + sample_processing_id,
                           df, FUN = sum, drop = FALSE)
  data$consensus_count[is.na(data$consensus_count)] <- 0

  meta <- unique(df[c("sequence_id", "length", "v_call", "j_call", "cdr3_aa")])
  data$length  <- meta$length [match(data$sequence_id, meta$sequence_id)]
  data$v_call  <- meta$v_call [match(data$sequence_id, meta$sequence_id)]
  data$j_call  <- meta$j_call [match(data$sequence_id, meta$sequence_id)]
  data$cdr3_aa <- meta$cdr3_aa[match(data$sequence_id, meta$sequence_id)]

  # Rename to match figure-5 column names
  data$seqid      <- data$sequence_id
  data$junction_aa <- data$cdr3_aa
  data$Read.count <- data$consensus_count

  new.data <- NULL
  for (z in unique(data$sample_processing_id)) {
    data.1    <- data[data$sample_processing_id == z, ]
    new.data.1 <- fig5_compute_hamming_by_length(data.1, threshold)
    new.data  <- rbind(new.data, new.data.1)
  }

  new.data
}

# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

test_that("fig5_update_read_counts matches package update_read_counts", {
  counts <- c(10, 5, 0, 3)
  sim    <- matrix(c(
    1.0, 0.9, 0.5, 0.3,
    0.9, 1.0, 0.6, 0.4,
    0.5, 0.6, 1.0, 0.8,
    0.3, 0.4, 0.8, 1.0
  ), nrow = 4)

  pkg_result <- update_counts(counts, sim, cutoff = 0.8)
  fig_result <- fig5_update_read_counts(counts, sim, threshold = 0.8)

  expect_equal(pkg_result, fig_result)
})

test_that("fig5 hamming similarity matches package hamming_similarity", {
  vj <- data.frame(
    cdr3_aa = c("ASGDKRGEQY", "ASGDSRNEQY", "ASGDTRGEQY"),
    junction_aa = c("ASGDKRGEQY", "ASGDSRNEQY", "ASGDTRGEQY"),
    length = 10,
    v_call = "TRBV13-2",
    j_call = "TRBJ2-7"
  )
  pkg_sim <- hamming_similarity(vj)

  dist.m <- stringdist::stringdistmatrix(vj$junction_aa, vj$junction_aa, method = "hamming")
  dist.m <- dist.m / 10
  fig_sim <- 1 - dist.m

  expect_equal(unname(pkg_sim), unname(fig_sim))
})

test_that("full pipeline: package and figure-5 logic return identical wrc values", {
  pkg_result <- adjust_counts(mouse_PBSvTCZ_data_minisubset,
                                      sim_method = "HAMMING",
                                      cutoff = 0.8)
  fig_result <- fig5_pipeline(mouse_PBSvTCZ_data_minisubset, threshold = 0.8)

  key_pkg <- paste(pkg_result$sequence_id, pkg_result$sample_processing_id)
  key_fig <- paste(fig_result$seqid,       fig_result$sample_processing_id)
  fig_wrc <- fig_result$wrc[match(key_pkg, key_fig)]

  expect_equal(pkg_result$wrc, fig_wrc, tolerance = 1e-9)
})
