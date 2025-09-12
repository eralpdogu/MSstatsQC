#' A function to convert mzQC files to MSstatsQC format
#'
#' @param mzqc data file to be converted
#' @return A data frame that can be used with MSstatsQC
#' @keywords mzqc, qcmetrics, input
#' @return A csv file that is converted from raw files
#' @export
#' @import jsonlite
#' @examples
#' \dontrun{library(RforProteomics)}
#' \dontrun{msfile <- getPXD000001mzXML()}
#' \dontrun{MSnbaseToMSstatsQC(msfile)}

mzQCToMSstatsQC  <-  function(mzQCfile) {
  
  mzqc_data <- fromJSON("mzqc-reanalyse.mzqc")
  
  MSstatsQCdata <- data.frame()
  
  for (i in 1:49) {
    df<-cbind(
      AcquiredTime = mzqc_data$mzQC$runQualities$metadata$inputFiles[[i]]$fileProperties[[1]]$value,
      Precursor = mzqc_data$mzQC$runQualities$qualityMetrics[[i]]$value[[2]]$Precursor,
      Annotations = mzqc_data$mzQC$runQualities$qualityMetrics[[i]]$value[[2]]$Annotations,
      RT = mzqc_data$mzQC$runQualities$qualityMetrics[[i]]$value[[2]]$BestRetentionTime,
      FWHM = mzqc_data$mzQC$runQualities$qualityMetrics[[i]]$value[[2]]$MaxFwhm,
      MinStartTime = mzqc_data$mzQC$runQualities$qualityMetrics[[i]]$value[[2]]$MinStartTime,
      MaxEndTime = mzqc_data$mzQC$runQualities$qualityMetrics[[i]]$value[[2]]$MaxEndTime,
      TotalArea = mzqc_data$mzQC$runQualities$qualityMetrics[[i]]$value[[2]]$TotalArea)
    MSstatsQCdata<-rbind(MSstatsQCdata,df)
  }
  
  MSstatsQCdata[,3:dim(merged_df)[2]] <- as.data.frame(lapply(merged_df[,3:dim(merged_df)[2]], 
                                                          function(x) {if (is.character(x)) {as.numeric(x)} else {x}}))
  MSstatsQCdata[,1]<-as.factor(merged_df[,1])
  MSstatsQCdata[,2]<-as.factor(merged_df[,2])
  return(MSstatsQCdata)
}


