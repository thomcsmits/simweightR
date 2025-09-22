#' Use TCRsimilift
#'
#' Convenience function to do all TCRsimilift data processing in one call.
#' The function checks the data to ensure format, prepares the data for processing,
#' runs net_update_data to return an updated dataframe with similarity-altered
#' counts, and offers the option to export results as .Rds and .csv files.
#'
#' @param df Input dataframe of AIRR formatted immunological data.
#' @param sim_method Either HAMMING or BLOSUM.
#' @param export_results Boolean, whether to automatically run the export_results function.
#' @param output_directory Name of output directory.
#' @param cutoff Minimum similarity score to consider two sequences neighbours. Between 0 and 1. Default 0.8 .
#'
#' @returns Returns dataframe with extra column of adjusted counts based on similarity.
#' @export
#'
#' @examples
#' results <- TCRsimilift_calculate(my_dataframe, sim_method="BLOSUM", export_results=TRUE, output_directory="my_outputs")
#'
TCRsimilift_calculate <- function(df,
                                  sim_method="HAMMING",
                                  export_results=FALSE,
                                  output_directory = "outputs",
                                  cutoff = 0.8) {

  datacheck(df)
  df2 <- dataprep(df)
  new.data <- net_update_data(df2, sim_method = sim_method, cutoff=cutoff)
  if (export_results) {
    export_outputs(new.data, output_directory = "outputs")
  }
  return(new.data)
}
