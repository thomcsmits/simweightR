#' Call compute_dist_matrix for each sample.
#'
#' @param data Dataframe of AIRR format immunological data.
#' @param sim_method HAMMING or BLOSUM.
#'
#' @returns Returns dataframe with adjusted counts.
#' @export
#'
net_update_data <- function(data, sim_method="HAMMING") {
  new.data <- c()
  for(sample_id in unique(data$sample_processing_id)){
    data.1 <- data[data$sample_processing_id == sample_id,]
    new.data.1 <- compute_dist_matrix(data.1, sample_id, sim_method = sim_method)
    new.data <- rbind(new.data, new.data.1)
  }
  return(new.data)
  # write.csv(new.data,
  # gzfile(paste0("data/updated-data-rc1-", sample_size, ".csv.gz")))
}
