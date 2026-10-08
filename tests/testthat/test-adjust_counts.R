test_that("weight_counts returns x with a weighted matrix", {
  x <- prepare_counts(mouse_pbs_srbc_mini)
  result <- weight_counts(x, threshold = 0.8)
  expect_true(is.list(result))
  expect_true(is.matrix(result$weighted))
  # counts and sequences are passed through unchanged
  expect_identical(result$counts, x$counts)
  expect_identical(result$sequences, x$sequences)
})

test_that("weight_counts output has one row per sequence and one column per sample", {
  result <- weight_counts(prepare_counts(mouse_pbs_srbc_mini), threshold = 0.8)
  # prepare_counts creates the full sequences x samples matrix (0 where a
  # sequence was not observed); weighted has the same rows and columns
  d <- mouse_pbs_srbc_mini
  n_seqs <- length(unique(paste(d$v_call, d$cdr3_aa, d$j_call)))
  n_samples <- length(unique(d$sample_processing_id))
  expect_equal(dim(result$weighted), c(n_seqs, n_samples))
  expect_identical(dimnames(result$weighted), dimnames(result$counts))
})

test_that("weight_counts weighted values are non-negative", {
  result <- weight_counts(prepare_counts(mouse_pbs_srbc_mini), threshold = 0.8)
  expect_true(all(result$weighted >= 0))
})

test_that("weight_counts works with BLOSUM method", {
  x <- prepare_counts(mouse_pbs_srbc_mini)
  result <- weight_counts(x, threshold = 0.9, method = "blosum")
  expect_true(is.matrix(result$weighted))
  expect_true(all(result$weighted >= 0))
  # same as passing blosum_similarity() as a custom function
  custom <- weight_counts(x, threshold = 0.9, method = function(s) blosum_similarity(s))
  expect_identical(result$weighted, custom$weighted)
  # and different from Hamming
  expect_false(isTRUE(all.equal(result$weighted,
                                weight_counts(x, threshold = 0.9)$weighted)))
})

test_that("weight_counts with stricter threshold does not raise errors", {
  result <- weight_counts(prepare_counts(mouse_pbs_srbc_mini), threshold = 0.99)
  expect_true(is.matrix(result$weighted))
})

test_that("weight_counts with permissive threshold does not raise errors", {
  result <- weight_counts(prepare_counts(mouse_pbs_srbc_mini), threshold = 0.5)
  expect_true(is.matrix(result$weighted))
})

test_that("weight_counts accepts a custom similarity function", {
  # identity matrix: each sequence is only similar to itself, so no
  # neighbours and weighted should equal counts
  identity_sim <- function(sequences) diag(length(sequences))
  result <- weight_counts(prepare_counts(mouse_pbs_srbc_mini), threshold = 0.8,
                          method = identity_sim)
  expect_true(is.matrix(result$weighted))
  expect_equal(result$weighted, result$counts)
})

test_that("weight_counts with threshold 1 leaves counts unchanged", {
  # sequences are unique within a group, so with threshold 1 no sequence has a
  # neighbour except itself, and weighted should equal counts
  result <- weight_counts(prepare_counts(mouse_pbs_srbc_mini), threshold = 1)
  expect_equal(result$weighted, result$counts)
})

test_that("weight_counts errors when custom function returns a non-matrix", {
  bad_sim <- function(sequences) as.data.frame(diag(length(sequences)))
  expect_error(weight_counts(prepare_counts(mouse_pbs_srbc_mini), threshold = 0.8,
                             method = bad_sim),
               "must return a matrix")
})

test_that("weight_counts errors when custom function returns wrong dimensions", {
  wrong_size_sim <- function(sequences) diag(length(sequences) + 1)
  expect_error(weight_counts(prepare_counts(mouse_pbs_srbc_mini), threshold = 0.8,
                             method = wrong_size_sim),
               "similarity matrix must be")
})

test_that("weight_counts errors when custom function returns values outside [0, 1]", {
  out_of_range_sim <- function(sequences) {
    m <- diag(length(sequences))
    m[1, 1] <- 2
    m
  }
  expect_error(weight_counts(prepare_counts(mouse_pbs_srbc_mini), threshold = 0.8,
                             method = out_of_range_sim),
               "between 0 and 1")
})

test_that("weight_counts errors on unrecognised method", {
  x <- prepare_counts(mouse_pbs_srbc_mini)
  expect_error(weight_counts(x, threshold = 0.8, method = "UNKNOWN"),
               "must be \"hamming\", \"blosum\" or a function")
  # names are case sensitive
  expect_error(weight_counts(x, threshold = 0.8, method = "HAMMING"),
               "must be \"hamming\", \"blosum\" or a function")
})

test_that("weight_counts matches a brute-force reference", {
  # For each sequence separately: find the sequences with the same V, J and
  # length, count mismatching letters, and apply the formula
  #   weighted_i = sum_j w_ij * count_j / sum_j w_ij,  w_ij = sim_ij if >= 0.8, else 0
  # This does not use split(), hamming_similarity() or matrix products, so it
  # also checks that every group's results end up in the right rows.
  x <- prepare_counts(mouse_pbs_srbc_mini)
  s <- x$sequences
  ref <- x$counts
  for (i in seq_len(nrow(s))) {
    same <- which(s$v_call == s$v_call[i] & s$j_call == s$j_call[i] &
                    nchar(s$cdr3_aa) == nchar(s$cdr3_aa[i]))
    letters_i <- strsplit(s$cdr3_aa[i], "")[[1]]
    sim <- sapply(s$cdr3_aa[same], function(b) {
      1 - sum(letters_i != strsplit(b, "")[[1]]) / length(letters_i)
    })
    w <- ifelse(sim >= 0.8 - 1e-10, sim, 0)
    ref[i, ] <- colSums(w * x$counts[same, , drop = FALSE]) / sum(w)
  }
  result <- weight_counts(x, threshold = 0.8)
  expect_equal(result$weighted, ref)
  # make sure the comparison is meaningful: weighting changes many rows
  expect_gt(sum(rowSums(abs(result$weighted - x$counts)) > 0), 10)
})

test_that("weight_counts does not depend on the order of the input rows", {
  # shuffle the rows of counts and sequences together: every sequence should
  # get the same weighted counts as before, found back by its id
  x <- prepare_counts(mouse_pbs_srbc_mini)
  result <- weight_counts(x, threshold = 0.8)

  set.seed(1)
  p <- sample(nrow(x$counts))
  shuffled <- list(counts = x$counts[p, ], sequences = x$sequences[p, ])
  result_shuffled <- weight_counts(shuffled, threshold = 0.8)

  expect_equal(result_shuffled$weighted[rownames(result$weighted), ], result$weighted)
})
