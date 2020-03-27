#' A function to convert MSnbase files to MSstatsQC format
#'
#' @param msfile data file to be converted
#' @return A data frame that can be used with MSstatsQC
#' @keywords MSnbase, qcmetrics, input
#' @return A csv file that is converted from raw files
#' @export
#' @importFrom MSnbase readMSData addIdentificationData fData rtime precursorIntensity
#' @import qcmetrics
#' @examples
#' \dontrun{library(RforProteomics)}
#' \dontrun{msfile <- getPXD000001mzXML()}
#' \dontrun{MSnbaseToMSstatsQC(msfile)}

MSnbaseToMSstatsQC  <-  function(msfile) {

  data <- readMSData(msfile, verbose = FALSE)

  if (!inherits(data, "MSnExp")) {
    stop("Only MSnSet class can be converted to input format for MSstats.")
  }
  qc <- QcMetric(name = "NULL")

  #Examples of metrics that can be monotired ###############################
  RetentionTime <- rtime(data)
  PrecursorIntensity <- precursorIntensity(data)
  ##########################################################################
  qcdata(qc, "RetentionTime") <- RetentionTime
  qcdata(qc, "PrecursorIntensity") <- PrecursorIntensity

  MSstatsQCdata <- c()
  MSstatsQCdata <- data.frame(setNames(lapply(ls(qc@qcdata), get, envir=qc@qcdata), ls(qc@qcdata)))
  MSstatsQCdata <- data.frame(AcquiredTime=seq_along(RetentionTime), Precursor=NA, Annotations=NA, MSstatsQCdata)

  ## if there are any missing variable name, warn it and stop
  check.name <- c("AcquiredTime", "Precursor", "Annotations", "RetentionTime", "PrecursorIntensity")

  diff.name <- setdiff(check.name, colnames(MSstatsQCdata))
  if (length(diff.name) > 0){
    stop(paste("Please check the variable name. The provided variable name", paste(diff.name, collapse=","), "is not present in the data set.", sep=" "))
  }
  return(MSstatsQCdata)
}

