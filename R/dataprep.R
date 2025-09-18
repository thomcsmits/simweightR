dataprep <- function(df) {
  data <- aggregate( consensus_count ~ sequence_id + sample_processing_id,
                     df, FUN=sum, drop=FALSE )
  data$consensus_count[is.na(data$consensus_count)] <- 0
  # data$treatment <- ifelse(grepl("PBS", data$sample_processing_id), "nt", "t") #REMOVE IN PKG

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
