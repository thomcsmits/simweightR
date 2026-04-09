#' Data preparation for simweightR
#'
#'
#' Preprocessing the data for downstream use involves:
#' * adds a column for cdr3 amino acid sequence lengths
#' * aggregates all counts for each sequence id and each sample.
#' * replaces NAs with 0.
#' * returns aggregate data with only necessary associated columns.
#'
#' See \link{adjust_counts} for a full example of the DGE workflow.
#'
#' @param df AIRR formatted dataframe of immunological data.
#'
#' @returns A dataframe
#'

dataprep <- function(df) {
  df$length <- nchar(df$cdr3_aa) #calculating cdr3 aa length

  # aggregate counts across different measurements wihin same sample
  data <- stats::aggregate( consensus_count ~ sequence_id + sample_processing_id,
                     df, FUN=sum, drop=FALSE )
  data$consensus_count[is.na(data$consensus_count)] <- 0

  uniqueIDVJ <- unique(df[c("sequence_id",
                            "length",
                            "v_call",
                            "j_call",
                            "cdr3_aa")])

  data$length <- uniqueIDVJ$length[match(data$sequence_id,
                                         uniqueIDVJ$sequence_id)]
  data$v_call <- uniqueIDVJ$v_call[match(data$sequence_id,
                                         uniqueIDVJ$sequence_id)]
  data$j_call <- uniqueIDVJ$j_call[match(data$sequence_id,
                                         uniqueIDVJ$sequence_id)]
  data$cdr3_aa <- uniqueIDVJ$cdr3_aa[match(data$sequence_id,
                                           uniqueIDVJ$sequence_id)]

  return(data)
}
