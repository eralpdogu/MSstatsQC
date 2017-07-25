#' A function to construct individual (X) and moving range (mR) control charts
#'
#' @param data comma-separated (.csv), metric file. It should contain a "Precursor" column and the metrics columns. It should also include "Annotations" for each observation.
#' @param peptide the name of precursor of interest.
#' @param L Lower bound of the guide set.
#' @param U Upper bound of the guide set.
#' @param metric the name of metric of interest.
#' @param normalization TRUE if metric is standardized and FALSE if not standardized.
#' @param ytitle the y-axis title of the plot.  Defaults to "Individual observations". The x-axis title is by default "Time : name of peptide"
#' @param type the type of the control chart. Two values can be assigned, "mean" or "variability". Default is "mean".
#' @param selectMean the mean of a metric. It is used when mean is known. It is NULL when mean is not known.  The default is NULL.
#' @param selectSD the standard deviation of a metric. It is used when standard deviation is known. It is NULL when mean is not known. The default is NULL.
#' @return A plot of individual values or moving ranges versus time per peptide and metric generated from \code{XmR.data.prepare} data frame.
#' @keywords XmR, control chart
#' @export
#' @import dplyr
#' @import RecordLinkage
#' @importFrom plotly plot_ly add_trace add_lines layout
#' @importFrom stats setNames sd
#' @examples
#' # First process the data to make sure it's ready to use
#' sampleData <- DataProcess(S9Site54)
#' head(sampleData)
#' # Find the name of the peptides
#' levels(sampleData$Precursor)
#' # Calculate X and mR statistics
#' XmRChart(data = sampleData, peptide = "VLVLDTDYK", metric = "BestRetentionTime")
#' XmRChart(data = sampleData, peptide = "VLVLDTDYK", metric = "BestRetentionTime",
#'          ytitle = "moving ranges", type = "variability")
#' XmRChart(data = sampleData, peptide = "VLVLDTDYK", metric = "BestRetentionTime",
#'          selectMean = 27.78, selectSD = 8.19)
#' XmRChart(data = sampleData, peptide = "DDGSWEVIEGYR", metric = "TotalArea")
#' XmRChart(data = sampleData, peptide = "DDGSWEVIEGYR", metric = "TotalArea",
#'          selectMean = 35097129, selectSD = 34132861)
#' XmRChart(data = sampleData, peptide = "TAAYVNAIEK", metric = "MaxFWHM")
#' XmRChart(data = sampleData, peptide = "LVNELTEFAK", metric = "Peak Assymetry")
################################################################################################################
XmRChart <- function(data = NULL, peptide, L = 1, U = 5, metric, normalization = FALSE,
                     ytitle = "Individual observations", type = "mean",
                     selectMean = NULL, selectSD = NULL) {
  #data <- input_checking(data)
  if(is.null(data))
    return()
  if(!is.data.frame(data)){
    stop(data)
  }
  metricData <- getMetricData(data, peptide, L, U, metric, normalization, selectMean, selectSD)
  precursor.data <- data[data$Precursor==peptide,]
  plot.data <- XmR.data.prepare(metricData, L, U, type, selectMean, selectSD)

  pal <- c("blue","red")
  pal <- setNames(pal,c("InRange","OutRange"))
  x <- list(
    title = paste("Time : ", peptide)
  )
  y <- list(
    title = ytitle
  )

  plot_ly(plot.data, x = ~QCno, y = ~t,showlegend = TRUE) %>%
    add_trace(x = ~QCno, y = ~t, color = ~InRangeOutRange, type="scatter",
              mode="markers", colors = pal , showlegend = TRUE) %>%
    add_lines(x = ~QCno, y = ~t, color = I("cornflowerblue"), showlegend = FALSE) %>%
    add_lines(y = ~LCL, color = I("red"), name = "LCL", showlegend = FALSE) %>%
    add_lines(y = ~UCL, color = I("red"), name = "UCL", showlegend = FALSE) %>%
    layout(xaxis = x,yaxis = y)
}
