#' A function to create heatmaps to compare performance with user defined performance criteria
#'
#' @param data Comma-separated (*.csv), QC file format. It should contain a Precursor
#'  column and the metrics columns.
#' @param method It is either "CUSUM" or "XmR"
#' @param peptideThresholdRed Is a threshold that marks percentage of peptides above it
#'  red on the heatmap. Defaults to 0.7
#' @param peptideThresholdYellow Is a threshold that marks percentage of peptides above
#'  it and below the peptideThresholdRed, yellow on the heatmap. Defaults to 0.5
#' @param L Lower bound of the giude set. Defaults to 1
#' @param U Upper bound of the guide set. Defaults to 5
#' @param listMean List of the means for the metrics. If you don't know the means leave
#'  it as NULL and they will be calculated automatically by using L and U. The default is NULL.
#' @param listSD List of the standard deviations for the metrics. If you don't know the
#'  standard deviations leave it as NULL and they will be calculated automatically by using L and U.
#'   The default is NULL.
#' @param type can take two values, "mean" or "dispersion". Defaults to "mean"
#' @param title the title of the plot. Defaults to "heatmap plot"
#' @return A heatmap to aggregate results per metric generated from \code{heatmap.DataFrame} data frame.
#' @keywords heatmap
#' @export
#' @import ggplot2
#' @import RecordLinkage
#' @import grid
#' @importFrom ggExtra removeGrid rotateTextX
#' @examples
#' # First process the data to make sure it's ready to use
#' sampleData <- DataProcess(S9Site54)
#' head(sampleData)
#' # Draw Decision maker plot
#' DecisionMap(data = sampleData, method = "CUSUM")
#' DecisionMap(data = sampleData, method = "CUSUM", type = "variability")
#' DecisionMap(data = sampleData, method = "XmR")
#' DecisionMap(data = sampleData, method = "XmR", type = "variability")
#' DecisionMap(data = sampleData, method = "CUSUM", type = "mean",
#'               listMean = list("BestRetentionTime" = 27.78,
#'                               "TotalArea" = 35097129,
#'                               "MaxFWHM" = 0.28,
#'                               "Peak Assymetry" = 0.98),
#'               listSD = list("BestRetentionTime" = 8.19,
#'                             "TotalArea" = 34132861,
#'                             "MaxFWHM" = 0.054,
#'                             "Peak Assymetry" = 0.002)
#'                  )

#########################################################################################################
DecisionMap <- function(data = NULL, method = "XmR",
                          peptideThresholdRed = 0.7, peptideThresholdYellow = 0.5,
                          L = 1, U = 5, type = "mean", title = "heatmap plot",
                          listMean = NULL, listSD = NULL) {

  if(is.null(data))
    return()
  if(!is.data.frame(data)){
    stop(data)
  }

  data.metrics <- c(find_custom_metrics(data))
  remove <- c("MinStartTime","MaxEndTime")
  data.metrics <- data.metrics[!data.metrics %in% remove]

  data <- heatmap.DataFrame(data, data.metrics,method,peptideThresholdRed,
                            peptideThresholdYellow, L, U, type,listMean, listSD)

  p <- ggplot(data,aes(data$time,data$metric, group = data$bin, fill = data$bin))
  p <- p + scale_fill_manual(values =
                               c("Pass" = "blue","Fail" = "red","Warning" = "yellow"))
  p <- p + geom_tile(colour="white",size=.1)
  p <- p + coord_equal()
  p <- p + removeGrid()
  p <- p + rotateTextX()
  p <- p + ggtitle(title,subtitle = "")
  p <- p + labs(x=NULL, y=NULL)
  p <- p +  theme(axis.text=element_text(size=12),legend.title = element_blank())
  p
}
