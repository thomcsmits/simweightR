#' Mouse TCR sequencing data, PBS vs SRBC
#'
#' TCR beta repertoires from the T cell zone (TCZ) of the spleen of 4 control
#' mice (injected with PBS) and 4 mice immunized with sheep red blood cells
#' (SRBC, sampled after 4 days). Only sequences with V gene TRBV2 are included,
#' so that all neighbours of each sequence are in the data. Sequences with a
#' read count below 2 were removed in the original data. The data is formatted
#' according to AIRR standards. See `data-raw/mouse_data.R` for how it was made.
#'
#' @format ## `mouse_pbs_srbc`
#' A data frame with 16,616 rows (one per clonotype per sample) and 6 columns:
#' \describe{
#'   \item{cdr3_aa}{CDR3 amino acid sequence.}
#'   \item{v_call}{V gene.}
#'   \item{j_call}{J gene.}
#'   \item{consensus_count}{Read count.}
#'   \item{sample_processing_id}{Sample identifier: mouse, zone and treatment.}
#'   \item{treatment}{PBS (control) or SRBC_4d (immunized).}
#' }
#' @source <https://www.biorxiv.org/content/10.1101/2025.03.28.645951v1>
"mouse_pbs_srbc"

#' Mouse TCR sequencing data, PBS vs SRBC - small example subset
#'
#' Subset of [mouse_pbs_srbc] with two V/J/length groups: CDR3 length 9 with
#' V gene TRBV2 and J gene TRBJ1-1 or TRBJ2-1. Small enough for quick examples
#' and tests, with several pairs of similar sequences in each group.
#'
#' @format ## `mouse_pbs_srbc_mini`
#' A data frame with 71 rows and the same 6 columns as [mouse_pbs_srbc].
#' @source <https://www.biorxiv.org/content/10.1101/2025.03.28.645951v1>
"mouse_pbs_srbc_mini"

#' BLOSUM62-derived distance between amino acids
#'
#' Distance between each pair of the 20 standard amino acids, derived from the
#' BLOSUM62 target frequencies q_ij (Henikoff & Henikoff, 1992). Each row of
#' q_ij is normalised to a probability vector (the probabilities of the amino
#' acids aligned to that amino acid); the distance between two amino acids is
#' the squared Euclidean distance between their vectors, scaled so that the
#' largest distance is 1, and rounded to 4 decimals. Used for the BLOSUM
#' similarity: 1 - (sum of per-position distances) / length. See
#' `data-raw/blosum/blosum62_distance.R`.
#'
#' @format ## `blosum62_distance`
#' A symmetric 20 x 20 numeric matrix with the one-letter amino acid codes as
#' row and column names (in BLOSUM order: A R N D C Q E G H I L K M F P S T W
#' Y V). The diagonal is 0; off-diagonal values range from 0.0456 (I-V) to 1.
#' @source Henikoff S, Henikoff JG (1992). Amino acid substitution matrices
#'   from protein blocks. PNAS 89(22):10915-10919.
#'   <https://doi.org/10.1073/pnas.89.22.10915>
"blosum62_distance"
