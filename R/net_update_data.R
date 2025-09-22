#' Call compute_dist_matrix for each sample.
#'
#' Creates a new.data empty dataframe. Then, for each sample, it updates the rows
#' of the input dataframe with adjusted counts, and appends them to new.data .
#' In the end, new.data is returned, consisting of all original input rows, but
#' with an added wrc column for similarity-based adjusted counts.
#'
#' @param data Dataframe of AIRR format immunological data.
#' @inheritParams TCRsimilift_calculate
#'
#' @returns Returns dataframe with adjusted counts.
#' @export
#'
net_update_data <- function(data, sim_method="HAMMING", cutoff=0.8) {
  new.data <- c()
  for(sample_id in unique(data$sample_processing_id)){
    data.1 <- data[data$sample_processing_id == sample_id,]
    new.data.1 <- compute_dist_matrix(data.1, sim_method = sim_method, cutoff=cutoff)
    new.data <- rbind(new.data, new.data.1)
  }
  return(new.data)
}
