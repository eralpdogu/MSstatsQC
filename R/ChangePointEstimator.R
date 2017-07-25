#' A function to identify the time of a change in the mean or variability of a metric
#'
#' @param data comma-separated (.csv), metric file. It should contain a "Precursor" column and the metrics columns. It should also include "Annotations" for each observation.
#' @param peptide the name of precursor of interest.
#' @param L Lower bound of the guide set.
#' @param U Upper bound of the guide set.
#' @param metric the name of metric of interest.
#' @param normalization TRUE metric is standardized and FALSE if not standardized.
#' @param ytitle the y-axis title of the plot.  Defaults to "Change Point Plot - mean". The x-axis title is by default "QCno-name of peptide"
#' @param type the type of the control chart. Two values can be assigned, "mean" or "variability". Default is "mean".
#' @param selectMean the mean of a metric. It is used when mean is known. It is NULL when mean is not known.  The default is NULL.
#' @param selectSD the standard deviation of a metric. It is used when standard deviation is known. It is NULL when mean is not known. The default is NULL.
#' @return A plot of likelihood statistics versus time per peptide and metric generated from \code{CP.data.prepare} data frame.
#' @keywords change point, control chart
#' @export
#' @importFrom plotly plot_ly add_markers add_lines layout
#' @import RecordLinkage
#' @import dplyr
#' @examples
#' # First process the data to make sure it's ready to use
#' sampleData <- DataProcess(S9Site54)
#' head(sampleData)
#' # Find the name of the peptides
#' levels(sampleData$Precursor)
#' # Calculate change point statistics
#' ChangePointEstimator(data = sampleData, peptide = "VLVLDTDYK", metric = "BestRetentionTime")
#' ChangePointEstimator(data = sampleData, peptide = "VLVLDTDYK", metric = "BestRetentionTime",
#'                      ytitle = "Change Point Plot - variability", type = "variability")
#' ChangePointEstimator(data = sampleData, peptide = "VLVLDTDYK", metric = "BestRetentionTime",
#'                      selectMean = 27.78, selectSD = 8.19)
#' ChangePointEstimator(data = sampleData, peptide = "DDGSWEVIEGYR", metric = "TotalArea")
#' ChangePointEstimator(data = sampleData, peptide = "DDGSWEVIEGYR", metric = "TotalArea",
#'                      selectMean = 35097129, selectSD = 34132861)
#' ChangePointEstimator(data = sampleData, peptide = "TAAYVNAIEK", metric = "MaxFWHM")
#' ChangePointEstimator(data = sampleData, peptide = "LVNELTEFAK", metric = "Peak Assymetry")

ChangePointEstimator <- function(data = NULL, peptide, L = 1, U = 5, metric, normalization = TRUE,
                                 ytitle = "Change Point Plot - mean", type = "mean", selectMean = NULL,
                                 selectSD = NULL) {
  if(is.null(data))
    return()
  if(!is.data.frame(data)){
    stop(data)
  }
  metricData <- getMetricData(data, peptide, L, U, metric, normalization, selectMean, selectSD)
  precursor.data <- data[data$Precursor==peptide,]
  ## Create variables
  plot.data <- CP.data.prepare(metricData, type)

  x <- list(
    title = paste("Time : ", peptide)
  )
  y <- list(
    title = ytitle
  )

  plot_ly(plot.data, x = ~QCno, y = ~Et,showlegend = FALSE)%>% #,text=precursor.data$Annotations)
    add_lines(x = ~tho.hat, color = I("red"))%>%
    add_lines(x = ~QCno, y = ~Et, color = I("cornflowerblue"))%>%
    add_markers(x = ~QCno, y = ~Et, color = I("blue"))%>%
    layout(xaxis = x,yaxis = y)
}
