packageBaseDir <- system.file(package = "MSstatsQC")
appDir <- paste0(packageBaseDir,"/shiny-examples/msstats-qc")
source(paste0(appDir,"/QCMetrics.R"))
source(paste0(appDir,"/ggradar.R"))
source(paste0(appDir,"/helper-functions.R"))

library(dplyr)
library(ggplot2)
library(scales)
library(grid)
library(tidyr)
library(viridis)
library(extrafont)

CUSUM.outrange.thld <- 5

#################################################################################################################
#INPUT : "prodata" is the data user uploads.
#        "precursorSelection" is the precursor that user selects in Data Import tab. it can be either one precursor(peptide) or it can be "all peptides"
#        "L" and "U" are lower and upper bound of guide set that user choose in Data Import tab.
#        "metric" is one of these metrics: COL.BEST.RET,COL.FWHM, COL.TOTAL.AREA,COL.PEAK.ASS or a metric that user defines in his data set
#        "plot.method" is one of XmR, CUSUM or CP methods
#        "normalization" is either TRUE or FALSE
#        "y.title1" and "y.title2" are titles of left and right plot which are Individual Value and Moving Range
# DESCRIPTION : Draws together the "Individual Value" and "Moving Range plots" (left and right plot) for each metric and method.
render.QC.chart <- function(prodata, precursorSelection, L, U, metric, plot.method, normalization, y.title1, y.title2,selectMean=NULL,selectSD=NULL, guidset_selected){
  validate(
    need(!is.null(prodata), "Please upload your data")
  )
  precursors <- levels(prodata$Precursor)
  plots <- list()
  annot_list <- list()
  if(precursorSelection == "all peptides") {
    results <- lapply(c(1:nlevels(prodata$Precursor)), function(j) {
      metricData <- getMetricData(prodata, precursors[j], L, U, metric = metric, normalization = normalization,selectMean,selectSD, guidset_selected)
      plots[[2*j-1]] <<- do.plot(prodata, metricData, precursors[j],L,U, plot.method, y.title1, type = 1,selectMean,selectSD, guidset_selected)
      plots[[2*j]] <<- do.plot(prodata, metricData, precursors[j],L,U, plot.method, y.title2, type = 2,selectMean,selectSD, guidset_selected)
      if(j==1) {
        plots[[2*j-1]] <<- plots[[2*j-1]] %>% layout(annotations = list(
          list(x = 0.5 , y = 1.05, text = "Mean", showarrow = F, xref='paper', yref='paper'),
          list(x = 0.38 , y = -0.2, text = precursors[j], showarrow = F, xref='paper', yref='paper'),
          list(x = 0.3, y = -0.2, text = "QCno - ", showarrow = F, xref = 'paper',yref='paper')
          #,list(x = -0.12, y = 1, text = "YLABEL", showarrow = F, xref = 'paper',yref='paper')
        ))
        plots[[2*j]] <<- plots[[2*j]] %>%
          layout(
            annotations = list(
              list(x = 0.5 , y = 1.05, text = "Dispersion", showarrow = F, xref='paper', yref='paper'),
              list(x = 0.65 , y = -0.2, text = precursors[j], showarrow = F, xref='paper', yref='paper'),
              list(x = 0.5, y = -0.2, text = "QCno - ", showarrow = F, xref = 'paper',yref='paper')
            ))
      } else {
        plots[[2*j-1]] <<- plots[[2*j-1]] %>% layout(annotations = list(
          list(x = 0.38 , y = -0.2, text = precursors[j], showarrow = F, xref='paper', yref='paper'),
          list(x = 0.3, y = -0.2, text = "QCno - ", showarrow = F, xref = 'paper',yref='paper')
          #,list(x = -0.13, y = 1, text = "YLABEL", showarrow = F, xref = 'paper',yref='paper')
        ))
        plots[[2*j]] <<- plots[[2*j]] %>%
          layout(
            annotations = list(
              list(x = 0.64 , y = -0.2, text = precursors[j], showarrow = F, xref='paper', yref='paper'),
              list(x = 0.5 , y = -0.2, text = "QCno - ", showarrow = F, xref='paper', yref='paper')
            ))
      }
    })

    do.call(subplot,c(plots,nrows=nlevels(prodata$Precursor))) %>%
      layout(autosize = F, width = 1400, height = nlevels(prodata$Precursor)*200)
  }

  else {
    metricData <- getMetricData(prodata, precursorSelection, L, U, metric = metric, normalization,selectMean,selectSD, guidset_selected)
    plot1 <- do.plot(prodata, metricData, precursorSelection,L,U, plot.method,  y.title1, type = 1,selectMean,selectSD, guidset_selected) %>%
      layout(
        annotations = list(
          list(x = 0.2 , y = 1.05, text = "Mean", showarrow = F, xref='paper', yref='paper'),
          list(x = 0.45 , y = -0.09, text = precursorSelection, showarrow = F, xref='paper', yref='paper'),
          list(x = 0.3 , y = -0.09, text = "QCno - ", showarrow = F, xref='paper', yref='paper')
          #,list(x = -0.13 , y = 1, text = "YLABEL", showarrow = F, xref='paper', yref='paper')
        )
      )
    plot2 <- do.plot(prodata, metricData, precursorSelection,L,U, plot.method,  y.title2, type = 2,selectMean,selectSD, guidset_selected) %>%
      layout(
        annotations = list(
          list(x = 0.8 , y = 1.05, text = "Dispersion", showarrow = F, xref='paper', yref='paper'),
          list(x = 0.8 , y = -0.09, text = precursorSelection, showarrow = F, xref='paper', yref='paper'),
          list(x = 0.5 , y = -0.09, text = "QCno - ", showarrow = F, xref='paper', yref='paper')

        )
      )

    subplot(plot1,plot2)
  }
}
#################################################################################################################
# INPUTS : "prodata" is the data user uploads.
#          "metricData" is the column of the data related to the metric we want. Forexample if we want retention time, it gives retention time column
#          "precursorSelection" is the precursor that user selects in Data Import tab. it can be either one precursor(peptide) or it can be "all peptides"
#          "L" and "U" are lower and upper bound of guide set that user choose in Data Import tab.
#          "plot.method" is one of XmR, CUSUM or CP methods
#          "y.title" is the title of the plot which is either Individual Value or Moving Range
#          "type" is either 1 or 2. one is "Individual Value" plot and other "Moving Range" plot
#DESCRIPTION : draw one plot (which is either Individual Value or Moving Range based on the type user chooses) for each metric and method
do.plot <- function(prodata, metricData, precursorSelection, L, U, plot.method,  y.title, type,selectMean,selectSD, guidset_selected) {
  if(plot.method=="CUSUM") {
    CUSUM.plot(prodata, metricData, precursorSelection, y.title, type)
  } else if(plot.method=="CP") {
    CP.plot(prodata, metricData, precursorSelection, y.title, type)
  } else if(plot.method=="XmR") {
    XmR.plot(prodata, metricData, precursorSelection, L, U, y.title, type,selectMean,selectSD, guidset_selected)
  }
}
#################################################################################################
# INPUTS : "prodata" is the data user uploads.
#          "metricData" is the column of the data related to the metric we want. Forexample if we want retention time, it gives retention time column
#          "precursorSelection" is the precursor that user selects in Data Import tab. it can be either one precursor(peptide) or it can be "all peptides"
#          "L" and "U" are lower and upper bound of guide set that user choose in Data Import tab.
#          "ytitle" is the title of the plot which is either Individual Value or Moving Range
#          "type" is either 1 or 2. one is "Individual Value" plot and other "Moving Range" plot
#DESCRIPTION: draws one CUSUM plot based on type for each given metric
CUSUM.plot <- function(prodata, metricData, precursorSelection,  ytitle, type) {
  plot.data <- CUSUM.data.prepare(prodata, metricData, precursorSelection, type)
  plot.data1 <- data.frame(
    QCno = rep(plot.data$QCno,2),
    CUSUMValue = c(plot.data$CUSUM.poz, plot.data$CUSUM.neg),
    Annotations = rep(plot.data$Annotations,2),
    outRangeInRange = c(as.character(plot.data$outRangeInRangePoz), as.character(plot.data$outRangeInRangeNeg))
  )
  #print(plot.data1)
  pal <- c("lightslateblue","red","blue","red")
  pal <- setNames(pal,c("InRangeCUSUM-","OutRangeCUSUM-","InRangeCUSUM+","OutRangeCUSUM+"))

  plot_ly(plot.data1, x = ~QCno, y = ~CUSUMValue,showlegend = FALSE, text = ~Annotations)%>%

    add_markers(x = ~QCno, y = ~CUSUMValue, color = ~outRangeInRange,
               type="scatter",mode="markers", colors = pal , showlegend = FALSE) %>%

    add_lines(y = CUSUM.outrange.thld, color = I("red"), name = "CUSUM thld", showlegend = FALSE) %>%
    add_lines(y = -CUSUM.outrange.thld, color = I("red"), name = "CUSUM thld", showlegend = FALSE) %>%
    add_lines(data = plot.data, x = ~QCno, y = ~CUSUM.poz, color = I("cornflowerblue"),name = "CUSUM+", showlegend = FALSE) %>%
    add_lines(data = plot.data, x = ~QCno, y = ~CUSUM.neg, color = I("lightskyblue"), name = "CUSUM-",showlegend = FALSE)

}
#########################################################################################################################
# INPUTS : "prodata" is the data user uploads.
#          "metricData" is the column of the data related to the metric we want. Forexample if we want retention time, it gives retention time column
#          "precursorSelection" is the precursor that user selects in Data Import tab. it can be either one precursor(peptide) or it can be "all peptides"
#          "ytitle" is the title of the plot which is either Individual Value or Moving Range
#          "type" is either 1 or 2. one is "Individual Value" plot and other "Moving Range" plot
#DESCRIPTION: draws one CP plot based on type for each given metric
CP.plot <- function(prodata, metricData, precursorSelection, ytitle, type) {
  precursor.data <- prodata[prodata$Precursor==precursorSelection,]
  if(type == 1) {Annotations = precursor.data$Annotations[-1]}
  if(type == 2) {Annotations = precursor.data$Annotations}
   plot.data <- CP.data.prepare(prodata, metricData, type)
   plot.data1 <- data.frame(
     QCno = plot.data$QCno,
     Et = plot.data$Et,
     tho.hat = plot.data$tho.hat,
     Annotations = Annotations
   )
    plot_ly(plot.data1, x = ~QCno, y = ~Et,showlegend = FALSE, text = ~Annotations)%>% #,text=precursor.data$Annotations)
      add_lines(x = ~tho.hat, color = I("red"))%>%
      add_lines(x = ~QCno, y = ~Et, color = I("cornflowerblue"))%>%
      add_markers(x = ~QCno, y = ~Et, color = I("blue"))
}
#########################################################################################################################
# INPUTS : "prodata" is the data user uploads.
#          "metricData" is the column of the data related to the metric we want. Forexample if we want retention time, it gives retention time column
#          "precursorSelection" is the precursor that user selects in Data Import tab. it can be either one precursor(peptide) or it can be "all peptides"
#          "L" and "U" are lower and upper bound of guide set that user choose in Data Import tab.
#          "ytitle" is the title of the plot which is either Individual Value or Moving Range
#          "type" is either 1 or 2. one is "Individual Value" plot and other "Moving Range" plot
#DESCRIPTION: draws one XmR plot based on type for each given metric
XmR.plot <- function(prodata, metricData, precursorSelection, L, U, ytitle, type,selectMean,selectSD, guidset_selected) {
  precursor.data <- prodata[prodata$Precursor==precursorSelection,]
  plot.data <- XmR.data.prepare(prodata, metricData, L, U, type,selectMean,selectSD, guidset_selected)
  plot.data1 <- data.frame(
    QCno = plot.data$QCno,
    t = plot.data$t,
    UCL = plot.data$UCL,
    LCL = plot.data$LCL,
    InRangeOutRange = plot.data$InRangeOutRange,
    Annotations = precursor.data$Annotations
  )

  pal <- c("blue","red")
  pal <- setNames(pal,c("InRange","OutRange"))

  plot_ly(plot.data1, x = ~QCno, y = ~t ,showlegend = FALSE, text = ~Annotations) %>%
    add_trace(x = ~QCno, y = ~t, color = ~InRangeOutRange, type="scatter",
              mode="markers", colors = pal , showlegend = FALSE) %>%
    add_lines(x = ~QCno, y = ~t, color = I("cornflowerblue"), showlegend = FALSE) %>%
    add_lines(y = ~LCL, color = I("red"), name = "LCL", showlegend = FALSE) %>%
    add_lines(y = ~UCL, color = I("red"), name = "UCL", showlegend = FALSE)

}
#################################################################################################################
XmR.Summary.plot <- function(prodata,data.metrics, L, U,listMean,listSD, guidset_selected) {
  dat <- XmR.Summary.DataFrame(prodata,data.metrics, L, U,listMean,listSD, guidset_selected)
  tho.hat.df <- get_CP_tho.hat(prodata, L, U, data.metrics,listMean,listSD, guidset_selected)
  gg <- ggplot(dat)
  gg <- gg + geom_hline(yintercept=0, alpha=0.5)
  gg <- gg + geom_smooth(method="loess",aes(x=dat$QCno, y=dat$pr.y,colour = group, group = group))
  gg <- gg + geom_point(data = tho.hat.df, aes(x = tho.hat.df$tho.hat, y = tho.hat.df$y, colour = "Change point"))
  gg <- gg + scale_color_manual(breaks = c("Metric mean increase",
                                           "Metric mean decrease",
                                           "Metric dispersion increase",
                                           "Metric dispersion decrease",
                                           "Change point"),
                                values = c("Metric mean increase" = "#E69F00",
                                           "Metric mean decrease" = "#56B4E9",
                                           "Metric dispersion increase" = "#009E73",
                                           "Metric dispersion decrease" = "#D55E00",
                                           "Change point" = "red"),
                                guide='legend')
  gg <- gg + guides(colour = guide_legend(override.aes = list(linetype=c(1,1,1,1,0),shape=c(NA,NA,NA,NA,16))))
  gg <- gg + facet_wrap(~metric,nrow = ceiling(length(data.metrics)/4))
  gg <- gg + annotate("text", x = 15, y = 1.3, label = "Mean")
  gg <- gg + annotate("text", x = 25, y = -1.3, label = "Dispersion")
  gg <- gg + scale_y_continuous(expand=c(0,0), limits = c(-1.4,1.4),breaks = c(1,0.5,0,-0.5,-1) ,labels = c(1,0.5,0,"0.5","1"))
  gg <- gg + labs(x = "QC No", y = "% of out of control \nprecursors")
  gg <- gg + ggtitle("Overall Summary \nXmR")
  theme_set(theme_gray(base_size = 15)) # this will change the size of all the texts in all ggplot functions
  gg <- gg + theme(plot.title = element_text(size=15, face="bold",margin = margin(10, 0, 10, 0)),
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
###############################################################################################
CUSUM.Summary.plot <- function(prodata, data.metrics, L, U,listMean,listSD, guidset_selected) {
   h <- 5
   dat <- CUSUM.Summary.DataFrame(prodata, data.metrics, L, U,listMean,listSD, guidset_selected)
   tho.hat.df <- get_CP_tho.hat(prodata, L, U, data.metrics,listMean,listSD, guidset_selected)

   gg <- ggplot(dat)
   gg <- gg + geom_hline(yintercept=0, alpha=0.5)
   gg <- gg + stat_smooth(method="loess", aes(x=dat$QCno, y=dat$pr.y, colour = group, group = group))
   gg <- gg + geom_point(data = tho.hat.df, aes(x = tho.hat.df$tho.hat, y = tho.hat.df$y, colour = "Change point"))
   #gg <- gg + geom_point(aes(x = c(16,6),y = c(1.1,-1.1), colour = "Change point"))
   gg <- gg + scale_color_manual(breaks = c("Metric mean increase",
                                            "Metric mean decrease",
                                            "Metric dispersion increase",
                                            "Metric dispersion decrease",
                                            "Change point"),
                                 values = c("Metric mean increase" = "#E69F00",
                                            "Metric mean decrease" = "#56B4E9",
                                            "Metric dispersion increase" = "#009E73",
                                            "Metric dispersion decrease" = "#D55E00",
                                            "Change point" = "red"),
                                 guide='legend')
   gg <- gg + guides(colour = guide_legend(override.aes = list(linetype=c(1,1,1,1,0),shape=c(NA,NA,NA,NA,16))))
   gg <- gg + facet_wrap(~metric,nrow = ceiling(length(data.metrics)/4))
   gg <- gg + annotate("text", x = 15, y = 1.3, label = "Mean")
   gg <- gg + annotate("text", x = 25, y = -1.3, label = "Dispersion")
   gg <- gg + scale_y_continuous(expand=c(0,0), limits = c(-1.4,1.4),
                                 breaks = c(1,0.5,0,-0.5,-1) ,labels = c(1,0.5,0,"0.5","1"))
   gg <- gg + ggtitle("Overall Summary \nCUSUM")

   gg <- gg + labs(x = "QC No", y = "% of out of control \nprecursors")
   gg <- gg + theme(plot.title = element_text(size=15, face="bold",margin = margin(10, 0, 10, 0)),
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
####################################################################
XmR.Radar.Plot <- function(prodata, data.metrics, L,U,listMean,listSD,guidset_selected) {

  dat <- XmR.Radar.Plot.DataFrame(prodata, data.metrics, L,U,listMean,listSD,guidset_selected)
  #write.csv(file="dataRadar.csv",dat)
  #print(dat)
  ggplot(dat, aes(y = OutRangeQCno, x = reorder(peptides,orderby),
                  group = group, colour = group, fill=group)) +
    coord_polar() +
    geom_point() +
    scale_fill_manual(breaks = c("Metric mean increase",
                                 "Metric mean decrease",
                                 "Metric dispersion increase",
                                 "Metric dispersion decrease"),
                      values = c("Metric mean increase" = "#E69F00",
                                 "Metric mean decrease" = "#56B4E9",
                                 "Metric dispersion increase" = "#009E73",
                                 "Metric dispersion decrease" = "#D55E00")) +
    scale_color_manual(breaks = c("Metric mean increase",
                                  "Metric mean decrease",
                                  "Metric dispersion increase",
                                  "Metric dispersion decrease"),
                       values = c("Metric mean increase" = "#E69F00",
                                  "Metric mean decrease" = "#56B4E9",
                                  "Metric dispersion increase" = "#009E73",
                                  "Metric dispersion decrease" = "#D55E00")) +
    facet_wrap(~metric,nrow = ceiling(length(data.metrics)/4)) +
    geom_polygon(alpha=0.5)+
    ggtitle("Precursor Level Summary \nXmR") +
    xlab("") +
    ylab("# of out of control \nQC samples") +
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

#################################################################################################################
CUSUM.Radar.Plot <- function(prodata, data.metrics, L,U,listMean,listSD,guidset_selected) {
  dat <- CUSUM.Radar.Plot.DataFrame(prodata, data.metrics, L,U,listMean,listSD,guidset_selected)

  ggplot(dat, aes(y = OutRangeQCno, x = reorder(peptides,orderby),
                  group = group, colour = group, fill = group)) +
    coord_polar() +
    geom_point() +
    scale_fill_manual(breaks = c("Metric mean increase",
                                 "Metric mean decrease",
                                 "Metric dispersion increase",
                                 "Metric dispersion decrease"),
                      values = c("Metric mean increase" = "#E69F00",
                                 "Metric mean decrease" = "#56B4E9",
                                 "Metric dispersion increase" = "#009E73",
                                 "Metric dispersion decrease" = "#D55E00")) +
    scale_color_manual(breaks = c("Metric mean increase",
                                  "Metric mean decrease",
                                  "Metric dispersion increase",
                                  "Metric dispersion decrease"),
                       values = c("Metric mean increase" = "#E69F00",
                                  "Metric mean decrease" = "#56B4E9",
                                  "Metric dispersion increase" = "#009E73",
                                  "Metric dispersion decrease" = "#D55E00")) +
    facet_wrap(~metric,nrow = ceiling(length(data.metrics)/4)) +
    geom_polygon(alpha=0.5)+
    ggtitle("Precursor Level Summary \nCUSUM") +
    xlab("") +
    ylab("# of out of control \nQC samples") +

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

#########################################################################################################################
metrics_box.plot <- function(prodata, data.metrics) {
  plots <- list()
  for(i in 1:length(data.metrics)) {
    metric <- data.metrics[i]
    precursor.data <- substring(reorder(prodata$Precursor,prodata[,metric]), first = 1, last = 3)
    plots[[i]] <- plot_ly(prodata, y = prodata[,metric], color = precursor.data, type = "box") %>%
      layout(
        annotations = list(
          list(x = 0.5 , y = 1, text = metric, showarrow = F, xref='paper', yref='paper')
        ),showlegend = FALSE)
      #layout(yaxis = list(title = metric),showlegend = FALSE)
  }
  height <- (length(data.metrics))*300
  p <- do.call(subplot,c(plots,nrows=length(plots))) %>%
    layout(autosize = F, width = 700, height = height)
  return(p)
}
#####################################################################################################
metrics_heat.map <- function(prodata,data.metrics, method,peptideThresholdRed,peptideThresholdYellow,
                             L, U, type, title,listMean, listSD, guidset_selected) {

  #color_palette <- colorRampPalette(c("green", "yellow", "red"))(3)
  data <- heatmap.DataFrame(prodata, data.metrics,method,peptideThresholdRed,peptideThresholdYellow,
                            L, U, type,listMean, listSD, guidset_selected)
  #print(data)
  p <- ggplot(data,aes(time,metric, group = flag, fill = flag))
  p <- p + scale_fill_manual(values=c("Acceptable" = "blue","Unacceptable" = "red","Poor" = "yellow")
                             )
  p <- p + geom_tile(colour="white",size=.1)
  p <- p + coord_equal()
  p <- p + removeGrid()
  p <- p + rotateTextX()
  p <- p + ggtitle(title,subtitle = "")
  p <- p + labs(x=NULL, y=NULL)
  p <- p +  theme(axis.text=element_text(size=12),legend.title = element_blank())
  p
}

