#' A function to create radar plot to aggregate results from X and mR charts or CUSUMm and CUSUMv charts.
#'
#' @param data omma-separated (.csv), metric file. It should contain a "Precursor" column and the metrics columns. It should also include "Annotations" for each observation.
#' @param L lower bound of the guide set.
#' @param U upper bound of the guide set.
#' @param method defines the method selected to construct control charts.
#' @param listMean list of the means for each metric. It is used when mean is known. It is NULL when mean is not known.  The default is NULL.
#' @param listSD list of the standard deviations for each metric. It is used when standard deviation is known. It is NULL when mean is not known. The default is NULL.
#' automatically by using L and U. The default is NULL.
#' @return A radar plot to aggregate results per metric generated from \code{XmR.Radar.Plot.DataFrame} data frame or \code{CUSUM.Radar.Plot.DataFrame} data frame.
#' @keywords XmR
#' @export
#' @import ggplot2
#' @import RecordLinkage
#' @import grid
#' @importFrom  stats reorder
#' @examples
#' # First process the data to make sure it's ready to use
#' sampleData <- DataProcess(S9Site54)
#' head(sampleData)
#' # Draw XmR radar plot
#' RadarPlot(data = sampleData)
#' RadarPlot(data = sampleData, method = "CUSUM")
#' RadarPlot(data = sampleData,
#'                 listMean = list("BestRetentionTime" = 27.78,
#'                                 "TotalArea" = 35097129,
#'                                 "MaxFWHM" = 0.28,
#'                                 "Peak Assymetry" = 0.98),
#'                 listSD = list("BestRetentionTime" = 8.19,
#'                               "TotalArea" = 34132861,
#'                               "MaxFWHM" = 0.054,
#'                               "Peak Assymetry" = 0.002)
#'                 )
###########################################################################################
RadarPlot <- function(data = NULL, L=1, U=5, method = "XmR", listMean=NULL, listSD=NULL) {

  if(method == "XmR") {
    if(is.null(data))
      return()
    data.metrics <- c(find_custom_metrics(data))
    remove <- c("MinStartTime","MaxEndTime")
    data.metrics <- data.metrics[!data.metrics %in% remove]
    dat <- XmR.Radar.Plot.DataFrame(data, data.metrics, L,U,listMean,listSD)
    OutRangeQCno <- dat$OutRangeQCno
    peptides <- dat$peptides
    orderby <- dat$orderby
    group <- dat$group
    ggplot(dat, aes(y = OutRangeQCno, x = reorder(peptides,orderby),
                    group = group, colour = group, fill=group)) +
      coord_polar() +
      geom_point() +
      scale_fill_manual(breaks = c("Mean increase",
                                   "Mean decrease",
                                   "Variability increase",
                                   "Variability decrease"),
                        values = c("Mean increase" = "#E69F00",
                                   "Mean decrease" = "#56B4E9",
                                   "Variability increase" = "#009E73",
                                   "Variability decrease" = "#D55E00")) +
      scale_color_manual(breaks = c("Mean increase",
                                    "Mean decrease",
                                    "Variability increase",
                                    "Variability decrease"),
                         values = c("Mean increase" = "#E69F00",
                                    "Mean decrease" = "#56B4E9",
                                    "Variability increase" = "#009E73",
                                    "Variability decrease" = "#D55E00")) +
      facet_wrap(~metric,nrow = ceiling(length(data.metrics)/4)) +
      geom_polygon(alpha=0.5)+
      ggtitle("Radar plot : X and mR") +
      xlab("") +
      ylab("# of out of control \nruns") +
      theme(
        axis.text.x = element_text(face="bold",size = rel(0.7)),
        axis.title.y=element_text(size=12),
        axis.text.y=element_text(size=12, hjust=0.5),
        plot.title = element_text(size=15, face="bold",margin = margin(10, 0, 10, 0)),
        legend.title=element_blank(),
        legend.text = element_text(size = 12),
        panel.grid.major = element_line(colour = "firebrick3",linetype = "dotted"),
        plot.margin = unit(c(1,3,1,1), "lines")
      )
  }

  else if(method == "CUSUM") {
    if(is.null(data))
      return()
    data.metrics <- c(find_custom_metrics(data))
    remove <- c("MinStartTime","MaxEndTime")
    data.metrics <- data.metrics[!data.metrics %in% remove]
    dat <- CUSUM.Radar.Plot.DataFrame(data, data.metrics, L,U,listMean,listSD)
    OutRangeQCno <- dat$OutRangeQCno
    peptides <- dat$peptides
    orderby <- dat$orderby
    group <- dat$group

    ggplot(dat, aes(y = OutRangeQCno, x = reorder(peptides,orderby),
                    group = group, colour = group, fill = group)) +
      coord_polar() +
      geom_point() +
      scale_fill_manual(breaks = c("Mean increase",
                                   "Mean decrease",
                                   "Variability increase",
                                   "Variability decrease"),
                        values = c("Mean increase" = "#E69F00",
                                   "Mean decrease" = "#56B4E9",
                                   "Variability increase" = "#009E73",
                                   "Variability decrease" = "#D55E00")) +
      scale_color_manual(breaks = c("Mean increase",
                                    "Mean decrease",
                                    "Variability increase",
                                    "Variability decrease"),
                         values = c("Mean increase" = "#E69F00",
                                    "Mean decrease" = "#56B4E9",
                                    "Variability increase" = "#009E73",
                                    "Variability decrease" = "#D55E00")) +
      facet_wrap(~metric,nrow = ceiling(length(data.metrics)/4)) +
      geom_polygon(alpha=0.5)+
      ggtitle("Radar plot : CUSUMm and CUSUMv") +
      xlab("") +
      ylab("# of out of control \nruns") +

      theme(
        axis.text.x = element_text(face="bold",size = rel(0.7)),
        axis.title.y=element_text(size=12),
        axis.text.y=element_text(size=12, hjust=0.5),
        plot.title = element_text(size=15, face="bold",margin = margin(10, 0, 10, 0)),
        legend.title=element_blank(),
        legend.text = element_text(size = 12),
        panel.grid.major = element_line(colour = "firebrick3",linetype = "dotted"),
        plot.margin = unit(c(1,3,1,1), "lines")
      )
  }

}

