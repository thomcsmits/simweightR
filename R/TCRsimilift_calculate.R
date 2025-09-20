TCRsimilift_calculate <- function(df, sim_method="HAMMING", export_results=FALSE) {
  datacheck(df)
  df2 <- dataprep(df)
  new.data <- net_update_data(df2, sim_method = sim_method)
  if (export_results) {
    smiley <- export_outputs(new.data)
  }
  return(new.data)
}
