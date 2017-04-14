#' A function to create cumulative sum charts for mean (CUSUMm) and cumulative sum charts for variability (CUSUMv) control charts
#'
#' @param data comma-separated (.csv), metric file. It should contain a "Precursor" column and the metrics columns. It should also include "Annotations" for each observation.
#' @param peptide the name of precursor of interest.
#' @param L Lower bound of the guide set.
#' @param U Upper bound of the guide set.
#' @param metric the name of metric of interest.
#' @param normalization TRUE if metric is standardized and FALSE if not standardized.
#' @param ytitle the y-axis title of the plot.  Defaults to "CUSUMm". The x-axis title is by default "Time : name of peptide"
#' @param type the type of the control chart. Two values can be assigned, "mean" or "dispersion". Default is "mean".
#' @param selectMean the mean of a metric. It is used when mean is known. It is NULL when mean is not known.  The default is NULL.
#' @param selectSD the standard deviation of a metric. It is used when standard deviation is known. It is NULL when mean is not known. The default is NULL.
#' @keywords cumulative Sum, control chart
#' @export
#' @import dplyr
#' @importFrom plotly plot_ly add_markers add_lines layout
#' @importFrom stats setNames
#' @import RecordLinkage
#' @examples
#' # First process the data to make sure it's ready to use
#' sampleData <- DataProcess(S9Site54)
#' head(sampleData)
#' # Find the name of the peptides
#' levels(sampleData$Precursor)
#' # Calculate CUSUM statistics
#' CUSUMPlot(data = sampleData, peptide = "VLVLDTDYK", metric = "BestRetentionTime")
#' CUSUMPlot(data = sampleData, peptide = "VLVLDTDYK", metric = "BestRetentionTime",
#'            ytitle = "CUSUMv", type = "dispersion")
#' CUSUMPlot(data = sampleData, peptide = "VLVLDTDYK", metric = "BestRetentionTime",
#'            selectMean = 27.78, selectSD = 8.19)
#' CUSUMPlot(data = sampleData, peptide = "DDGSWEVIEGYR", metric = "TotalArea")
#' CUSUMPlot(data = sampleData, peptide = "DDGSWEVIEGYR", metric = "TotalArea",
#'            selectMean = 35097129, selectSD = 34132861)
#' CUSUMPlot(data = sampleData, peptide = "TAAYVNAIEK", metric = "MaxFWHM")
#' CUSUMPlot(data = sampleData, peptide = "LVNELTEFAK", metric = "Peak Assymetry")

#################################################################################################
CUSUMChart<- function(data = NULL, peptide, L = 1, U = 5, metric, normalization = TRUE,
                      ytitle = "CUSUMm", type = "mean", selectMean = NULL, selectSD = NULL) {
  if(is.null(data))
    return()
  #data <- input_checking(data)
  if(!is.data.frame(data)){
    stop(data)
  }
  CUSUM.outrange.thld <- 5
  metricData <- getMetricData(data, peptide, L, U, metric, normalization, selectMean, selectSD)
  plot.data <- CUSUM.data.prepare(data, metricData, peptide, type)
  plot.data1 <- data.frame(
    QCno = rep(plot.data$QCno,2),
    CUSUMValue = c(plot.data$CUSUM.poz, plot.data$CUSUM.neg),
    Annotations = rep(plot.data$Annotations,2),
    outRangeInRange = c(as.character(plot.data$outRangeInRangePoz),
                        as.character(plot.data$outRangeInRangeNeg))
  )

  pal <- c("lightslateblue","red","blue","red")
  pal <- setNames(pal,c("InRangeCUSUM-","OutRangeCUSUM-","InRangeCUSUM+","OutRangeCUSUM+"))
  x <- list(
    title =  paste("Time : ", peptide),
    range = c(0, max(plot.data$QCno))
  )
  y <- list(
    title = ytitle
  )

  plot_ly(plot.data1, x = ~QCno, y = ~CUSUM.poz,showlegend = FALSE)%>%
    add_markers(x = ~QCno, y = ~CUSUMValue, color = ~outRangeInRange,
                type="scatter",mode="markers", colors = pal , showlegend = TRUE) %>%

    add_lines(y = CUSUM.outrange.thld, color = I("red"), name = "CUSUM thld", showlegend = FALSE) %>%
    add_lines(y = -CUSUM.outrange.thld, color = I("red"), name = "CUSUM thld", showlegend = FALSE) %>%
    add_lines(data = plot.data, x = ~QCno, y = ~CUSUM.poz, color = I("cornflowerblue"),
              name = "CUSUM+", showlegend = FALSE) %>%
    add_lines(data = plot.data, x = ~QCno, y = ~CUSUM.neg, color = I("lightskyblue"),
              name = "CUSUM-",showlegend = FALSE)%>%
    layout(xaxis = x,yaxis = y)
}
