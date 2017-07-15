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
#' @import ggplot2
#' @import RecordLinkage
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
    
  if(is.null(data))
    return()
  if(!is.data.frame(data)){
    stop(data)
  }
  data.metrics <- c(find_custom_metrics(data))
  remove <- c("MinStartTime","MaxEndTime")
  data.metrics <- data.metrics[!data.metrics %in% remove]
  dat <- XmR.Summary.DataFrame(data,data.metrics, L, U, listMean, listSD)
  tho.hat.df <- get_CP_tho.hat(data, L, U, data.metrics, listMean, listSD)

  gg <- ggplot(dat)
  gg <- gg + geom_hline(yintercept=0, alpha=0.5)
  gg <- gg + geom_smooth(method="loess",aes(x=dat$QCno, y=dat$pr.y,colour = dat$group,
                                            group = dat$group))
  gg <- gg + geom_point(data = tho.hat.df, aes(x = tho.hat.df$tho.hat, y = tho.hat.df$y,
                                               colour = "Change point"))
  gg <- gg + scale_color_manual(breaks = c("Mean increase",
                                           "Mean decrease",
                                           "Variability increase",
                                           "Variability decrease",
                                           "Change point"),
                                values = c("Mean increase" = "#E69F00",
                                           "Mean decrease" = "#56B4E9",
                                           "Variability increase" = "#009E73",
                                           "Variability decrease" = "#D55E00",
                                           "Change point" = "red"),
                                guide='legend')
  gg <- gg + guides(colour = guide_legend(override.aes = list(linetype=c(1,1,1,1,0),shape=c(NA,NA,NA,NA,16))))
  gg <- gg + facet_wrap(~metric,nrow = ceiling(length(data.metrics)/4))
  gg <- gg + annotate("text", x = 15, y = 1.3, label = "Mean")
  gg <- gg + annotate("text", x = 25, y = -1.3, label = "Variability")
  gg <- gg + scale_y_continuous(expand=c(0,0), limits = c(-1.4,1.4),
                                breaks = c(1,0.5,0,-0.5,-1) ,labels = c(1,0.5,0,"0.5","1"))
  gg <- gg + labs(x = "Time", y = "% of out of control \nruns")
  gg <- gg + ggtitle("River plots : X and mR")
  theme_set(theme_gray(base_size = 15)) # this will change the size of all
  #the texts in all ggplot functions
  gg <- gg + theme(plot.title = element_text(size=15, face="bold",
                                             margin = margin(10, 0, 10, 0)),
                   axis.text.x=element_text(size=12, vjust=0.5),
                   axis.text.y=element_text(size=12, hjust=0.5),
                   axis.title.y=element_text(size=12),
                   axis.title.x=element_text(size=12),
                   legend.text = element_text(size = 12),
                   legend.title=element_blank(),
                   plot.margin = unit(c(1,3,1,1), "lines")
  )
  gg

  }
  
  else if(method == "CUSUM") {
    
    if(is.null(data))
      return()
    if(!is.data.frame(data)){
      stop(data)
    }
    h <- 5
    data.metrics <- c(find_custom_metrics(data))
    remove <- c("MinStartTime","MaxEndTime")
    data.metrics <- data.metrics[!data.metrics %in% remove]
    
    dat <- CUSUM.Summary.DataFrame(data, data.metrics, L, U,listMean,listSD, decisionInterval=5)
    tho.hat.df <- get_CP_tho.hat(data, L, U, data.metrics,listMean,listSD)
    
    gg <- ggplot(dat)
    gg <- gg + geom_hline(yintercept=0, alpha=0.5)
    gg <- gg + stat_smooth(method="loess", aes(x=dat$QCno, y=dat$pr.y,
                                               colour = dat$group, group = dat$group))
    gg <- gg + geom_point(data = tho.hat.df,
                          aes(x = tho.hat.df$tho.hat, y = tho.hat.df$y,
                              colour = "Change point"))
    gg <- gg + scale_color_manual(breaks = c("Mean increase",
                                             "Mean decrease",
                                             "Variability increase",
                                             "Variability decrease",
                                             "Change point"),
                                  values = c("Mean increase" = "#E69F00",
                                             "Mean decrease" = "#56B4E9",
                                             "Variability increase" = "#009E73",
                                             "Variability decrease" = "#D55E00",
                                             "Change point" = "red"),
                                  guide='legend')
    gg <- gg + guides(colour =
                        guide_legend(
                          override.aes = list(linetype=c(1,1,1,1,0),shape=c(NA,NA,NA,NA,16))))
    gg <- gg + facet_wrap(~metric,nrow = ceiling(length(data.metrics)/4))
    gg <- gg + annotate("text", x = 15, y = 1.3, label = "Mean")
    gg <- gg + annotate("text", x = 25, y = -1.3, label = "Variability")
    gg <- gg + scale_y_continuous(expand=c(0,0), limits = c(-1.4,1.4),
                                  breaks = c(1,0.5,0,-0.5,-1) ,
                                  labels = c(1,0.5,0,"0.5","1"))
    gg <- gg + ggtitle("River plots : CUSUMm and CUSUMv")
    
    gg <- gg + labs(x = "Time", y = "% of out of control \npeptides")
    gg <- gg + theme(plot.title = element_text(size=15, face="bold",
                                               margin = margin(10, 0, 10, 0)),
                     axis.text.x=element_text(size=12, vjust=0.5),
                     axis.text.y=element_text(size=12, hjust=0.5),
                     axis.title.y=element_text(size=12),
                     axis.title.x=element_text(size=12),
                     legend.text = element_text(size = 12),
                     legend.title=element_blank(),
                     plot.margin = unit(c(1,3,1,1), "lines")
    )
    
    gt <- ggplot_gtable(ggplot_build(gg))
    gt$layout$clip[gt$layout$name == "panel"] <- "off"
    grid.draw(gt)
    gg
    
  }
}

