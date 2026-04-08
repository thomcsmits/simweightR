#' Use TCRsimilift
#'
#' Convenience function to do all TCRsimilift data processing in one call.
#' The function checks the data to ensure format, prepares the data for processing,
#' runs weight_counts to return an updated dataframe with similarity-altered
#' counts, and offers the option to export results as .Rds or .csv files.
#'
#' @param df Input dataframe of AIRR formatted immunological data.
#' @param sim_method Either HAMMING or BLOSUM.
#' @param export_results Boolean, whether to automatically run the export_results function.
#' @param output_directory Name of output directory.
#' @param cutoff Minimum similarity score to consider two sequences neighbours. Between 0 and 1. Default 0.8 .
#' @inheritParams export_outputs
#'
#'
#' @returns Returns dataframe with extra column of adjusted counts based on similarity.
#' @export
#'
#' @examples
#' results <- adjust_counts(mouse_PBSvTCZ_data_minisubset, sim_method="HAMMING", cutoff = 0.77, export_results=FALSE)
#'
#'#'@examples
#'\dontrun{
#' # This is an example of DGE incorporating our data preprocessing.
#' results <- adjust_counts(mouse_PBSvTCZ_data)
#' count_matrix <- as_counts_matrix(results, doFilter = TRUE)
#' # Differential gene expression analysis using Wilcoxon test.
#' # The functions DGEList(), calcNormFactors() and cpm() need the library edgeR.
#' class(count_matrix) <- "numeric"
#' sample_size <- 2
#' condition <- factor(c( rep("treatment", sample_size), rep("non-treatment", sample_size)))
#' y <- edgeR::DGEList( counts=count_matrix, group=condition )
#' y <- edgeR::calcNormFactors(y,method="TMM")
#' count_norm=edgeR::cpm(y)
#' count_norm<-as.data.frame(count_norm)
#' pvalues <- sapply(1:nrow(count_norm),function(i){
#'   data<-cbind.data.frame(gene=as.numeric(t(count_norm[i,])),condition)
#'   p=wilcox.test(gene~condition, data)$p.value
#'   return(p)
#'   })
#' pvalues <- data.frame(pvalues)
#' rownames(pvalues) <- rownames(count_norm)
#' pvalues$fdr=p.adjust(pvalues$pvalues,method = "fdr")
#' conditionsLevel<-levels(condition)
#' dataCon1=count_norm[,c(which(condition==conditionsLevel[1]))]
#' dataCon2=count_norm[,c(which(condition==conditionsLevel[2]))]
#' dataCon2 <- dataCon2[rownames(dataCon1), ]
#' pvalues <- pvalues[rownames(dataCon1), ]
#' foldChanges=log2((rowMeans(dataCon2) + 1) / (rowMeans(dataCon1) + 1))
#' outRst<-data.frame(log2foldChange=foldChanges, pValues=pvalues$pvalues, FDR=pvalues$fdr)
#' rownames(outRst)=rownames(count_norm)
#' outRst=na.omit(outRst)
#' fdrThres=1
#' tbl <- outRst[outRst$FDR<=fdrThres,]
#'}
#'
adjust_counts <- function(df,
                          sim_method="HAMMING",
                          cutoff = 0.8,
                          export_results=FALSE,
                          output_directory = "outputs",
                          csv_output = FALSE) {

  datacheck(df)
  df2 <- dataprep(df)
  new.data <- weight_counts(df2, sim_method = sim_method, cutoff = cutoff)
  if (export_results) {
    export_outputs(new.data, output_directory = output_directory, csv_output = csv_output)
  }
  else {
    return(new.data)
  }
}
