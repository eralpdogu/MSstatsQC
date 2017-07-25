#' A function to create river plot to aggregate results from X and mR charts or CUSUMm and CUSUMv charts.
#'
#' @param data omma-separated (.csv), metric file. It should contain a "Precursor" column and the metrics columns. It should also include "Annotations" for each observation.
#' @param L lower bound of the guide set.
#' @param U upper bound of the guide set.
#' @param listMean list of the means for each metric. It is used when mean is known. It is NULL when mean is not known.  The default is NULL.
#' @param listSD list of the standard deviations for each metric. It is used when standard deviation is known. It is NULL when mean is not known. The default is NULL.
#' @return A river plot to aggregate results per metric generated from \code{XmR.Summary.DataFrame} data frame or \code{CUSUM.Summary.DataFrame} data frame.
#' @keywords XmR
#' @export
#' @import RecordLinkage
#' @import ggplot2
#' @import grid
#' @examples
#' # First process the data to make sure it's ready to use
#' sampleData <- DataProcess(S9Site54)
#' head(sampleData)
#' # Draw XmR summary plot
#' RiverPlot(data = sampleData)
#' RiverPlot(data = sampleData, L=1, U=20, method = "XmR",
#'                 listMean = list("BestRetentionTime" = 27.78,
#'                                 "TotalArea" = 35097129,
#'                                  "MaxFWHM" = 0.28,
#'                                  "PeakAssymetry" = 0.98),
#'                 listSD = list("BestRetentionTime" = 8.19,
#'                               "TotalArea" = 34132861,
#'                               "MaxFWHM" = 0.054,
#'                               "PeakAssymetry" = 0.002)
#'                 )

RiverPlot <- function(data = NULL, L=1, U=5, method = "XmR", listMean=NULL, listSD=NULL) {

  if(method == "XmR") {
    gg <- SummaryPlot(data , L , U , method = "XmR",
                      listMean=NULL, listSD=NULL)
    gg <- gg + ggtitle("Overall Summary \nXmR")
    gg
  }

  else if(method == "CUSUM") {
    gg <- SummaryPlot(data , L , U , method = "CUSUM",
                      listMean=NULL, listSD=NULL)
    gg <- gg + ggtitle("Overall Summary \nCUSUM")
    gg
  }

}



