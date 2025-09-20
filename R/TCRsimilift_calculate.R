TCRsimilift_calculate <- function(df,
                                  sim_method="HAMMING",
                                  export_results=FALSE,
                                  output_directory = "outputs") {

  datacheck(df)
  df2 <- dataprep(df)
  new.data <- net_update_data(df2, sim_method = sim_method)
  if (export_results) {
    export_outputs(new.data, output_directory = "outputs")
  }
  return(new.data)
}
