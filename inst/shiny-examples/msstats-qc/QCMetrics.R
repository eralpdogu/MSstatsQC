#COL.BEST.RET <- "Retention Time"
#COL.FWHM <- "Full Width at Half Maximum"
#COL.TOTAL.AREA <- "Total Peak Area"
#COL.PEAK.ASS <- "Peak Assymetry"
#############################################################################################
#INPUTS : "prodata" is the data user uploads.
#         "precursorSelection" is the precursor that user selects in Data Import tab. it can be either one precursor(peptide) or it can be "all peptides"
#         "L" and "U" are lower and upper bound of guide set that user choose in Data Import tab.
#         "metric" is one of these metrics: COL.BEST.RET,COL.FWHM, COL.TOTAL.AREA,COL.PEAK.ASS or a metric that user defines in his data set
#         "normalization" is either TRUE or FALSE
#Description of function : it gets the metric column only for the precursor chosen and either return the column as it is or normalize it and then return it
getMetricData <- function(prodata, precursorSelection, L, U, metric, normalization,selectMean,selectSD, guidset_selected) {
  precursor.data<-prodata[prodata$Precursor==precursorSelection,] #"Precursor" is one of the columns in data that shows the name of peptides. This gets only the part of data related to the selected peptide
  metricData <- 0

  if(is.null(metric)){
    return(NULL)
  }

  metricData = precursor.data[,metric]
  if(normalization == TRUE) {
    if(guidset_selected) {
      mu=mean(metricData[L:U]) # in-control process mean
      sd=sd(metricData[L:U]) # in-control process variance

    }else {
      mu = selectMean
      sd = selectSD

    }

    if(sd == 0) {sd <- 0.0001}
    metricData=scale(metricData[1:length(metricData)],mu,sd) # transformation for N(0,1) )
    return(metricData)
  } else if(normalization == FALSE){
    return(metricData)
  }

}
#########################################################################################################
find_custom_metrics <- function(prodata) {

    if(is.null(prodata))
      return (NULL)
    prodata <- prodata[, which(colnames(prodata)=="Annotations"):ncol(prodata),drop = FALSE]
    nums <- sapply(prodata, is.numeric)
    other.metrics <- colnames(prodata[,nums])[1:ifelse(length(colnames(prodata[,nums]))<11,
                                                       length(colnames(prodata[,nums])),
                                                       10)
                                              ] # limiting custom metrics up to 10 metrics and not more
    if(any(is.na(other.metrics))) {
      return(c())
    }
    return(other.metrics)
}
################################################################
#z = metricData
#INPUT : "prodata" is the data user uploads.
#        "metricData" is the column of the data related to the metric we want. Forexample if we want retention time, it gives retention time column
#        "precursorSelection" is the precursor that user selects in Data Import tab. it can be either one precursor(peptide) or it can be "all peptides"
#        "L" and "U" are lower and upper bound of guide set that user choose in Data Import tab.
#        "type" is either 1 or 2. one is "Individual Value" plot and other "Moving Range" plot
#DESCRIPTION : returns a data frame for CUSUM that contains all the information needed to plot CUSUM
CUSUM.data.prepare <- function(prodata, metricData, precursorSelection, type) {

  k=0.5
  CUSUM.outrange.thld <- 5
  outRangeInRangePoz = rep(0,length(metricData))
  outRangeInRangeNeg = rep(0,length(metricData))
  #precursor.data only gets the data for the selected precursor
  precursor.data <- prodata[prodata$Precursor==precursorSelection,] #"Precursor" is one of the columns in data that shows the name of peptides

  v <- numeric(length(metricData))

  Cpoz <- numeric(length(metricData))
  Cneg <- numeric(length(metricData))

  for(i in 2:length(metricData)) {
    Cpoz[i]=max(0,(metricData[i]-(k)+Cpoz[i-1]))
    Cneg[i]=max(0,((-k)-metricData[i]+Cneg[i-1]))
  }

  if(type==2) {
    for(i in 2:length(metricData)) {
      v[i]=(sqrt(abs(metricData[i]))-0.822)/0.349
    }
    for(i in 2:length(metricData)) {
      Cpoz[i]=max(0,(v[i]-(k)+Cpoz[i-1]))
      Cneg[i]=max(0,((-k)-v[i]+Cneg[i-1]))
    }
  }

  QCno = 1:length(metricData)
  for(i in 1:length(metricData)) {
    if(Cpoz[i] >= CUSUM.outrange.thld || Cpoz[i] <= -CUSUM.outrange.thld)
      outRangeInRangePoz[i] <- "OutRangeCUSUM+"
    else
      outRangeInRangePoz[i] <- "InRangeCUSUM+"
  }
  for(i in 1:length(metricData)) {
    if(-Cneg[i] >= CUSUM.outrange.thld || -Cneg[i] <= -CUSUM.outrange.thld)
      outRangeInRangeNeg[i] <- "OutRangeCUSUM-"
    else
      outRangeInRangeNeg[i] <- "InRangeCUSUM-"
  }
  plot.data =
    data.frame(QCno = QCno
               ,CUSUM.poz = Cpoz
               ,CUSUM.neg = -Cneg
               ,Annotations=precursor.data$Annotations
               ,outRangeInRangePoz = outRangeInRangePoz
               ,outRangeInRangeNeg = outRangeInRangeNeg
               )

  return(plot.data)
}
###################################################################################################
#z = metricData
#INPUT : "prodata" is the data user uploads.
#        "metricData" is the column of the data related to the metric we want. Forexample if we want retention time, it gives retention time column
#        "type" is either 1 or 2. one is "Individual Value" plot and other "Moving Range" plot
#DESCRIPTION : returns a data frame for CP that contains all the information needed to plot Change Point
CP.data.prepare <- function(prodata, metricData, type) {
  length_metricData <- length(metricData)
  Et <-  numeric(length_metricData-1) # this is Ct in type 1, and Dt in type 2.
  SS<- numeric(length_metricData-1)
  SST<- numeric(length_metricData-1)
  tho.hat <- 0

  if(type == 1) {
    ## Change point analysis for mean (Single step change model)
    for(i in 1:(length_metricData - 1)) {
      Et[i]=(length_metricData-i)*(((1/(length_metricData-i))*sum(metricData[(i+1):length_metricData]))-0)^2 #change point function
    }
    QCno=1:(length_metricData - 1)
  } else if(type == 2) {
    ## Change point analysis for variance (Single step change model)
    for(i in 1:length_metricData) {
      SS[i]=metricData[i]^2
    }
    for(i in 1:length_metricData) {
      SST[i]=sum(SS[i:length_metricData])
      Et[i]=((SST[i]/2)-((length_metricData-i+1)/2)*log(SST[i]/(length_metricData-i+1))-(length_metricData-i+1)/2) #change point function
    }
    QCno=1:length_metricData
  }
  tho.hat = which(Et==max(Et)) # change point estimate
  plot.data <- data.frame(QCno,Et,tho.hat)

  return(plot.data)
}
###################################################################################################
#INPUT : "prodata" is the data user uploads.
#        "L" and "U" are lower and upper bound of guide set that user choose in Data Import tab.
#        "data.metrics" is all the available metrics. It is defined in server.R
get_CP_tho.hat <- function(prodata, L, U, data.metrics,listMean,listSD, guidset_selected) {
  tho.hat <- data.frame(tho.hat = c(), metric = c(), group = c(), y=c())
  #precursors <- levels(reorder(prodata$Precursor,prodata[,COL.BEST.RET]))
  precursors <- levels(prodata$Precursor)
  for(metric in data.metrics) {
    for (j in 1:nlevels(prodata$Precursor)) {

      metricData <- getMetricData(prodata, precursors[j], L, U, metric = metric, normalization = TRUE,selectMean = listMean[[metric]],selectSD = listSD[[metric]], guidset_selected)
      mix <- rbind(
        data.frame(tho.hat = CP.data.prepare(prodata, metricData, type = 1)$tho.hat[1], metric = metric, group = "Individual Value", y=1.1),
        data.frame(tho.hat = CP.data.prepare(prodata, metricData, type = 2)$tho.hat[1], metric = metric, group = "Moving Range", y=-1.1)
      )

      tho.hat <- rbind(tho.hat, mix)

    }
  }

  return(tho.hat)
}
###################################################################################################
#z = metricData
#INPUT : "prodata" is the data user uploads.
#        "metricData" is the column of the data related to the metric we want. Forexample if we want retention time, it gives retention time column
#        "L" and "U" are lower and upper bound of guide set that user choose in Data Import tab.
#        "type" is either 1 or 2. one is "Individual Value" plot and other "Moving Range" plot
#DESCRIPTION : returns a data frame for XmR that contains all the information needed to plot XmR
XmR.data.prepare <- function(prodata, metricData, L,U, type,selectMean,selectSD, guidset_selected) {
  t <- numeric(length(metricData)-1)
  UCL <- 0
  LCL <- 0
  InRangeOutRange <- rep(0,length(metricData))

  for(i in 2:length(metricData)) {
    t[i]=abs(metricData[i]-metricData[i-1]) # Compute moving range of metricData
  }

  QCno=1:length(metricData)

  if(type == 1) {
    if(guidset_selected) {
      UCL=mean(metricData[L:U])+2.66*sd(t[L:U])
      LCL=mean(metricData[L:U])-2.66*sd(t[L:U])
    }else {
      UCL = selectMean + 2.66 * selectSD
      LCL = selectMean - 2.66 * selectSD
    }
    t <- metricData

  } else if(type == 2) {
    ## Calculate MR chart statistics and limits
    if(guidset_selected) {
      UCL=3.267*sd(t[1:L-U])
    }else{
      UCL = 3.267 * selectSD
    }
    LCL=0
  }
  for(i in 1:length(metricData)) {
    if(t[i] > LCL && t[i] < UCL)
      InRangeOutRange[i] <- "InRange"
    else
      InRangeOutRange[i] <- "OutRange"
  }
  plot.data=data.frame(QCno,metricData,t,UCL,LCL,InRangeOutRange)

  return(plot.data)
}
############################################################################################
#INPUTS : "prodata" is the data user uploads.
#         "metric" is one of these metrics: COL.BEST.RET,COL.FWHM, COL.TOTAL.AREA,COL.PEAK.ASS or a metric that user defines in his data set
#         "L" and "U" are lower and upper bound of guide set that user choose in Data Import tab.
#        "type" is either 1 or 2. one is "Individual Value" plot and other "Moving Range" plot
#DESCRIPTION : returns a data frame that is used in CUSUM.Summary.DataFrame function below to use for summary plot of CUSUM in Summary Tab of shiny app
CUSUM.Summary.prepare <- function(prodata, metric, L, U,type,selectMean,selectSD, guidset_selected) {
  h <- 5

  QCno <- 1:nrow(prodata)
  y.poz <- rep(0,nrow(prodata))
  y.neg <- rep(0,nrow(prodata))
  counter <- rep(0,nrow(prodata))

  #precursors <- levels(reorder(prodata$Precursor,prodata[,COL.BEST.RET]))
  precursors <- levels(prodata$Precursor)

  for(j in 1:length(precursors)) {
    metricData <- getMetricData(prodata, precursors[j], L, U, metric = metric, normalization = T,selectMean,selectSD, guidset_selected)

    counter[1:length(metricData)] <- counter[1:length(metricData)]+1
    plot.data <- CUSUM.data.prepare(prodata, metricData, precursors[j], type)

    sub.poz <- plot.data[plot.data$CUSUM.poz >= h | plot.data$CUSUM.poz <= -h, ]
    sub.neg <- plot.data[plot.data$CUSUM.neg >= h | plot.data$CUSUM.neg <= -h, ]

    y.poz[sub.poz$QCno] <- y.poz[sub.poz$QCno] + 1
    y.neg[sub.neg$QCno] <- y.neg[sub.neg$QCno] + 1
  }
  max_QCno <- max(which(counter!=0))
  pr.y.poz = y.poz[1:max_QCno]/counter[1:max_QCno]
  pr.y.neg = y.neg[1:max_QCno]/counter[1:max_QCno]

  plot.data <- data.frame(QCno = rep(1:max_QCno,2),
                          pr.y = c(pr.y.poz, pr.y.neg),
                          group = ifelse(rep(type==1,2*max_QCno),
                                         c(rep("Metric mean increase",max_QCno),
                                           rep("Metric mean decrease",max_QCno)),
                                         c(rep("Metric dispersion increase",max_QCno),
                                           rep("Metric dispersion decrease",max_QCno))),
                          metric = rep(metric,max_QCno*2)
     )
  return(plot.data)
}
############################################################################################
CUSUM.Summary.DataFrame <- function(prodata, data.metrics, L, U,listMean,listSD, guidset_selected) {
  dat <- data.frame(QCno = c(),
                    pr.y = c(),
                    group = c(),
                    metric = c())
  for (metric in data.metrics) {
    data.1   <- CUSUM.Summary.prepare(prodata, metric = metric, L, U,type = 1,selectMean = listMean[[metric]],selectSD = listSD[[metric]], guidset_selected)
    data.2   <- CUSUM.Summary.prepare(prodata, metric = metric, L, U,type = 2,selectMean = listMean[[metric]],selectSD = listSD[[metric]], guidset_selected)
    data.2$pr.y <- -(data.2$pr.y)
    dat <- rbind(dat,data.1,data.2)
  }
  return(dat)
}
############################################################################################
#DESCRIPTION : for each metric returns a data frame of QCno, probability of out of control peptide for dispersion or mean plot
XmR.Summary.prepare <- function(prodata, metric, L, U,type,selectMean,selectSD, guidset_selected) {
  QCno    <- 1:nrow(prodata)
  y.poz <- rep(0,nrow(prodata))
  y.neg <- rep(0,nrow(prodata))
  counter <- rep(0,nrow(prodata))

  #precursors <- levels(reorder(prodata$Precursor,prodata[,COL.BEST.RET]))
  precursors <- levels(prodata$Precursor)

  for(j in 1:length(precursors)) {
    metricData <- getMetricData(prodata, precursors[j], L = L, U = U, metric = metric, normalization = T,selectMean,selectSD, guidset_selected)
    counter[1:length(metricData)] <- counter[1:length(metricData)]+1
    plot.data <- XmR.data.prepare(prodata, metricData , L , U , type, selectMean,selectSD, guidset_selected)

    sub.poz <- plot.data[plot.data$t >= plot.data$UCL, ]
    sub.neg <- plot.data[plot.data$t <= plot.data$LCL, ]

    y.poz[sub.poz$QCno] <- y.poz[sub.poz$QCno] + 1
    y.neg[sub.neg$QCno] <- y.neg[sub.neg$QCno] + 1
  }
  max_QCno <- max(which(counter!=0))
  pr.y.poz = y.poz[1:max_QCno]/counter[1:max_QCno]
  pr.y.neg = y.neg[1:max_QCno]/counter[1:max_QCno]

  plot.data <- data.frame(QCno = rep(1:max_QCno,2),
                          pr.y = c(pr.y.poz, pr.y.neg),
                          group = ifelse(rep(type==1,2*max_QCno),
                                         c(rep("Metric mean increase",max_QCno),
                                           rep("Metric mean decrease",max_QCno)),
                                         c(rep("Metric dispersion increase",max_QCno),
                                           rep("Metric dispersion decrease",max_QCno))),
                          metric = rep(metric,max_QCno*2))

  return(plot.data)
}
###########################################################################################
#DESCRIPTION : returns a data frame for all the metrics, probability of out of range peptides, wheter it is for metric mean or metric dispersion and i is an increase or decrease
XmR.Summary.DataFrame <- function(prodata, data.metrics, L, U,listMean, listSD, guidset_selected) {
  dat <- data.frame(QCno = c(),
                    pr.y = c(),
                    group = c(),
                    metric = c())
  for (metric in data.metrics) {
    data.1   <- XmR.Summary.prepare(prodata, metric = metric, L, U,type = 1,selectMean = listMean[[metric]],selectSD = listSD[[metric]], guidset_selected)
    data.2   <- XmR.Summary.prepare(prodata, metric = metric, L, U,type = 2,selectMean = listMean[[metric]],selectSD = listSD[[metric]], guidset_selected)
    data.2$pr.y <- -(data.2$pr.y)
    dat <- rbind(dat, data.1, data.2)
  }
  return(dat)
}
############################################################################################
heatmap.DataFrame <- function(prodata, data.metrics,method,peptideThresholdRed,peptideThresholdYellow,
                              L, U, type,listMean, listSD, guidset_selected) {
  time <- c()
  val <- c()
  met <- c()
  flag <- c()

  for (metric in data.metrics) {
    df <- Decision.DataFrame.prepare(prodata, metric, method, peptideThresholdRed,peptideThresholdYellow,
                                     L, U,type,selectMean=listMean[[metric]],selectSD=listSD[[metric]], guidset_selected)
    time_df <- as.character(df$AcquiredTime)
    val_df <- df$pr.y
    met_df <- rep(metric,length(val_df))
    flag_df <- df$flag
    time <- c(time,time_df)
    val <- c(val,val_df)
    met <- c(met,met_df)
    flag <- c(flag,flag_df)
  }

  dataFrame <- data.frame(time = time,
                          value = val,
                          metric = met,
                          flag = flag
                         )

  return(dataFrame)
}
############################################################################################
# Compute.QCno.OutOfRangePeptide.XmR <- function(prodata,L,U,metric,type, XmR.type,selectMean,selectSD, guidset_selected) {
#   #precursors <- levels(reorder(prodata$Precursor,prodata[,COL.BEST.RET]))
#   precursors <- levels(prodata$Precursor)
#   QCno.out.range <- c()
#
#   for(j in 1:length(precursors)) {
#     metricData <- getMetricData(prodata, precursors[j], L = L, U = U, metric = metric, normalization = T,selectMean,selectSD, guidset_selected)
#     plot.data <- XmR.data.prepare(prodata, metricData , L = L, U = U, type ,selectMean,selectSD, guidset_selected)
#     if(XmR.type == "poz")
#       QCno.out.range <- c(QCno.out.range,length(plot.data[plot.data$t >= plot.data$UCL, ]$QCno))
#     else
#       QCno.out.range <- c(QCno.out.range,length(plot.data[plot.data$t <= plot.data$LCL, ]$QCno))
#   }
#   return(QCno.out.range)
# }
#############################################################################################
# Compute.QCno.OutOfRangePeptide.CUSUM <- function(prodata,L,U,metric,type, CUSUM.type,selectMean,selectSD, guidset_selected) {
#   h <- 5
#   #precursors <- levels(reorder(prodata$Precursor,prodata[,COL.BEST.RET]))
#   precursors <- levels(prodata$Precursor)
#   QCno.out.range <- c()
#
#   for(j in 1:length(precursors)) {
#     metricData <- getMetricData(prodata, precursors[j], L, U, metric = metric, normalization = T,selectMean,selectSD, guidset_selected)
#     plot.data <- CUSUM.data.prepare(prodata, metricData, precursors[j], type)
#     if(CUSUM.type == "poz")
#       QCno.out.range <- c(QCno.out.range,length(plot.data[plot.data$CUSUM.poz >= h | plot.data$CUSUM.poz <= -h, ]$QCno))
#     else
#       QCno.out.range <- c(QCno.out.range,length(plot.data[plot.data$CUSUM.neg >= h | plot.data$CUSUM.neg <= -h, ]$QCno))
#   }
#   return(QCno.out.range)
# }
###############################################################################################################
XmR.Radar.Plot.prepare <- function(prodata,L,U, metric, type,group,XmR_type,selectMean,selectSD, guidset_selected) {
  precursors <- levels(prodata$Precursor)
  precursors2 <- substring(precursors, first = 1, last = 3)
  QCno.length <- c()
  QCno.out.range.poz <- c()
  QCno.out.range.neg <- c()
  for(j in 1:length(precursors)) {
    metricData <- getMetricData(prodata, precursors[j], L = L, U = U, metric = metric, normalization = T,selectMean,selectSD, guidset_selected)
    QCno.length <- c(QCno.length,length(metricData))
    plot.data <- XmR.data.prepare(prodata, metricData , L = L, U = U, type ,selectMean,selectSD, guidset_selected)
    QCno.out.range.poz <- c(QCno.out.range.poz,length(plot.data[plot.data$t >= plot.data$UCL, ]$QCno))
    QCno.out.range.neg <- c(QCno.out.range.neg,length(plot.data[plot.data$t <= plot.data$LCL, ]$QCno))
  }
  #pr <- Compute.QCno.OutOfRangePeptide.XmR(prodata,L,U,metric = metric,type = type, XmR.type,selectMean,selectSD, guidset_selected)
  if(XmR_type == "poz") {
    dat <- data.frame(peptides = precursors2,
                          OutRangeQCno  = QCno.out.range.poz,
                          group         = rep(group,length(precursors)),
                          orderby       = seq(1:length(precursors)),
                          metric        = rep(metric, length(precursors)),
                          tool          = rep("XmR",length(precursors)),
                          probability   = QCno.out.range.poz/QCno.length
    )
  } else {
    dat <- data.frame(peptides = precursors2,
                          OutRangeQCno  = QCno.out.range.neg,
                          group         = rep(group,length(precursors)),
                          orderby       = seq(1:length(precursors)),
                          metric        = rep(metric, length(precursors)),
                          tool          = rep("XmR",length(precursors)),
                          probability   = QCno.out.range.neg/QCno.length
    )
  }

   return(dat)
  #return(list(dat.poz,dat.neg))
}
################################################################################################
XmR.Radar.Plot.DataFrame <- function(prodata, data.metrics, L,U,listMean,listSD,guidset_selected) {
  dat <- data.frame(peptides = c(), OutRangeQCno = c(), group = c(),
                    orderby = c(), metric = c(), tool = c(),
                    probability   = c()
                    )
  for (metric in data.metrics) {
    data.1 <- XmR.Radar.Plot.prepare(prodata,L,U,metric = metric, type = 1,group = "Metric mean increase",XmR_type = "poz",
                                     selectMean = listMean[[metric]],selectSD = listSD[[metric]],guidset_selected)
    data.2 <- XmR.Radar.Plot.prepare(prodata,L,U,metric = metric, type = 1,group = "Metric mean decrease",XmR_type = "neg",
                                     selectMean = listMean[[metric]],selectSD = listSD[[metric]],guidset_selected)


    data.3 <- XmR.Radar.Plot.prepare(prodata,L,U,metric = metric, type = 2,group = "Metric dispersion increase",XmR_type = "poz",
                                     selectMean = listMean[[metric]],selectSD = listSD[[metric]],guidset_selected)
    data.4 <- XmR.Radar.Plot.prepare(prodata,L,U,metric = metric, type = 2,group = "Metric dispersion decrease",XmR_type = "neg",
                                     selectMean = listMean[[metric]],selectSD = listSD[[metric]],guidset_selected)

    dat <- rbind(dat, data.1, data.2, data.3, data.4)
  }

  return(dat)
}
#################################################################################################################
CUSUM.Radar.Plot.prepare <- function(prodata,L,U, metric,type,group, CUSUM.type,selectMean,selectSD, guidset_selected) {
  #precursors <- levels(reorder(prodata$Precursor,prodata[,COL.BEST.RET]))
  h <- 5
  precursors <- levels(prodata$Precursor)
  precursors2 <- substring(precursors, first = 1, last = 3)
  QCno.length <- c()
  QCno.out.range <- c()
  for(j in 1:length(precursors)) {
    metricData <- getMetricData(prodata, precursors[j], L = L, U = U, metric = metric, normalization = T,selectMean,selectSD, guidset_selected)
    QCno.length <- c(QCno.length,length(metricData))
    plot.data <- CUSUM.data.prepare(prodata, metricData, precursors[j], type)
    if(CUSUM.type == "poz")
      QCno.out.range <- c(QCno.out.range,length(plot.data[plot.data$CUSUM.poz >= h | plot.data$CUSUM.poz <= -h, ]$QCno))
    else
      QCno.out.range <- c(QCno.out.range,length(plot.data[plot.data$CUSUM.neg >= h | plot.data$CUSUM.neg <= -h, ]$QCno))
  }
  #pr <- Compute.QCno.OutOfRangePeptide.CUSUM(prodata,L,U,metric = metric,type = type, CUSUM.type,selectMean,selectSD, guidset_selected)
  dat <- data.frame(peptides = precursors2,
                    OutRangeQCno  = QCno.out.range,
                    group         = rep(group,length(precursors)),
                    orderby       = seq(1:length(precursors)),
                    metric        = rep(metric, length(precursors)),
                    tool          = rep("XmR",length(precursors)),
                    probability   = QCno.out.range/QCno.length
  )
  return(dat)
}
#################################################################################################
CUSUM.Radar.Plot.DataFrame <- function(prodata, data.metrics, L,U,listMean,listSD,guidset_selected) {
  dat <- data.frame(peptides = c(), OutRangeQCno = c(), group = c(),
                    orderby = c(), metric = c(), tool = c(),
                    probability   = c()
  )
  for (metric in data.metrics) {
   data.1 <- CUSUM.Radar.Plot.prepare(prodata,L,U, metric = metric, type = 1, group = "Metric mean increase", CUSUM.type = "poz",selectMean = listMean[[metric]],selectSD = listSD[[metric]],guidset_selected)
   data.2 <- CUSUM.Radar.Plot.prepare(prodata,L,U, metric = metric, type = 1, group = "Metric mean decrease", CUSUM.type = "neg",selectMean = listMean[[metric]],selectSD = listSD[[metric]],guidset_selected)
   data.3 <- CUSUM.Radar.Plot.prepare(prodata,L,U, metric = metric, type = 2, group = "Metric dispersion increase", CUSUM.type = "poz",selectMean = listMean[[metric]],selectSD = listSD[[metric]],guidset_selected)
   data.4 <- CUSUM.Radar.Plot.prepare(prodata,L,U, metric = metric, type = 2, group = "Metric dispersion decrease", CUSUM.type = "neg",selectMean = listMean[[metric]],selectSD = listSD[[metric]],guidset_selected)
   dat <- rbind(dat, data.1, data.2, data.3, data.4)
  }

  return(dat)
}
#######################################################################################################
Decision.DataFrame.prepare <- function(prodata, metric, method, peptideThresholdRed,peptideThresholdYellow,
                                       L, U,type,selectMean,selectSD, guidset_selected) {
  h <- 5
  AcquiredTime <- prodata$AcquiredTime
  QCno    <- 1:nrow(prodata)
  y <- rep(0,nrow(prodata))
  counter <- rep(0,nrow(prodata))

  #precursors <- levels(reorder(prodata$Precursor,prodata[,COL.BEST.RET]))
  precursors <- levels(prodata$Precursor)
  if(method == "XmR") {
    for(precursor in precursors) {
      metricData <- getMetricData(prodata, precursor, L = L, U = U, metric = metric, normalization = T,selectMean,selectSD, guidset_selected)
      counter[1:length(metricData)] <- counter[1:length(metricData)]+1
      plot.data <- XmR.data.prepare(prodata, metricData , L , U , type,selectMean,selectSD, guidset_selected)
      sub <- plot.data[plot.data$InRangeOutRange == "OutRange",]
      #sub1 <- plot.data[(plot.data$t >= plot.data$UCL), ]
      #sub2 <- plot.data[plot.data$t <= plot.data$LCL, ]
      #sub <- rbind(sub1,sub2)
      y[sub$QCno] <- y[sub$QCno] + 1
    }
  } else if(method == "CUSUM") {
    for(precursor in precursors) {
      metricData <- getMetricData(prodata, precursor, L = L, U = U, metric = metric, normalization = T,selectMean,selectSD, guidset_selected)
      counter[1:length(metricData)] <- counter[1:length(metricData)]+1
      plot.data <- CUSUM.data.prepare(prodata, metricData, precursor, type)
      sub <- plot.data[(plot.data$CUSUM.poz >= h | plot.data$CUSUM.poz <= -h) | (plot.data$CUSUM.neg >= h | plot.data$CUSUM.neg <= -h), ]
      #sub2 <- plot.data[plot.data$CUSUM.neg >= h | plot.data$CUSUM.neg <= -h, ]
      #sub <- rbind(sub1,sub2)
      y[sub$QCno] <- y[sub$QCno] + 1
    }
  }
  max_QCno <- max(which(counter!=0))

  pr.y = y[1:max_QCno]/counter[1:max_QCno]


    plot.data <- data.frame(AcquiredTime = AcquiredTime[1:max_QCno],
                            QCno = rep(1:max_QCno,1),
                            pr.y = pr.y,
                            group = ifelse(rep(type==1,max_QCno),
                                           rep("Metric mean",max_QCno),
                                           rep("Metric dispersion",max_QCno)
                            ),
                            metric = rep(metric,max_QCno),

                            flag = rep(0,max_QCno)

                            )



  for (i in 1:max_QCno) {
    if(plot.data$pr.y[i] > peptideThresholdRed){
      plot.data$flag[i] <- "Unacceptable"
    }
    else if(plot.data$pr.y[i] > peptideThresholdYellow){
      plot.data$flag[i] <- "Poor"
    }
    else {
      plot.data$flag[i] <- "Acceptable"
    }
  }

    if(type == 2) {

      return(plot.data[-1,])
    }
    return(plot.data)



}
#######################################################################################################
 number.Of.Out.Of.Range.Metrics <- function(prodata,data.metrics,method, peptideThresholdRed,
                                            peptideThresholdYellow, L, U, type,listMean,listSD,
                                            guidset_selected) {
  #change this part
  metricCounterAboveRed = 0
  metricCounterAboveYellowBelowRed = 0
  precursors <- levels(prodata$Precursor)
  #### check prodata$Precursor vs levels(preata$Precursor) in for
  for (metric in data.metrics) {

  QCno    <- 1:nrow(prodata)
  y <- rep(0,nrow(prodata))
  counter <- rep(0,nrow(prodata))
  for(precursor in precursors) {
    metricData <- getMetricData(prodata, precursor, L = L, U = U, metric = metric, normalization = T,listMean[[metric]],listSD[[metric]], guidset_selected)
    counter[1:length(metricData)] <- counter[1:length(metricData)]+1
    #if(method == "XmR") {
      plot.data <- XmR.data.prepare(prodata, metricData , L , U , type,selectMean = listMean[[metric]],selectSD = listSD[[metric]], guidset_selected)
    #}

    sub <- plot.data[plot.data$InRangeOutRange == "OutRange",]
    y[sub$QCno] <- y[sub$QCno] + 1
  }
  max_QCno <- max(which(counter!=0))
  pr.y = y[1:max_QCno]/counter[1:max_QCno]

  if(type == 2) {
    pr.y <- pr.y[-1]
  }

  aboveYellow <- which(pr.y > peptideThresholdYellow)
  aboveYellowBelowRed <- which(pr.y > peptideThresholdYellow & pr.y <= peptideThresholdRed)

  if(length(which(pr.y > peptideThresholdRed)) > 0) {
    metricCounterAboveRed = metricCounterAboveRed + 1
  }
  if(length(aboveYellowBelowRed) > 0) {
    metricCounterAboveYellowBelowRed = metricCounterAboveYellowBelowRed + 1

  }
  }
  return(c(metricCounterAboveRed,metricCounterAboveYellowBelowRed))
}
####################################################################################################
decisionRule_warning_message_XmR <- function(prodata, data.metrics, method, peptideThresholdRed, peptideThresholdYellow,metricThresholdRed,metricThresholdYellow, L,U, type, listMean,listSD , guidset_selected) {


  a1 <- number.Of.Out.Of.Range.Metrics(prodata,data.metrics,method = "XmR", peptideThresholdRed,peptideThresholdYellow,
                                      L, U, type = 1,listMean,listSD, guidset_selected)
  XmRCounterAboveRed1 <- a1[1]
  XmRCounterAboveYellow1 <- a1[2]

  a2 <- number.Of.Out.Of.Range.Metrics(prodata,data.metrics,method = "XmR", peptideThresholdRed,peptideThresholdYellow,
                                       L, U, type = 2,listMean,listSD, guidset_selected)
  XmRCounterAboveRed2 <- a2[1]
  XmRCounterAboveYellow2 <- a2[2]

  if(XmRCounterAboveRed1 > metricThresholdRed && XmRCounterAboveRed2 > metricThresholdRed) return({paste("XmR RED FLAG: System performance is UNACCEPTABLE")})
  if(XmRCounterAboveRed1 > metricThresholdRed && XmRCounterAboveRed2 <= metricThresholdRed) return({paste("XmR RED FLAG: System performance is UNACCEPTABLE")})
  if(XmRCounterAboveRed1 <= metricThresholdRed && XmRCounterAboveRed2 > metricThresholdRed) return({paste("XmR RED FLAG: System performance is UNACCEPTABLE")})
  if(XmRCounterAboveYellow1 > metricThresholdYellow && XmRCounterAboveYellow2 > metricThresholdYellow) return({"XmR Yellow FLAG: System performance is POOR"})
  if(XmRCounterAboveYellow1 <= metricThresholdYellow && XmRCounterAboveYellow2 > metricThresholdYellow) return({"XmR Yellow FLAG: System performance is POOR"})
  if(XmRCounterAboveYellow1 > metricThresholdYellow && XmRCounterAboveYellow2 <= metricThresholdYellow) return({"XmR Yellow FLAG: System performance is POOR"})
  return({"XmR BLUE FLAG: System performance is acceptable"})
}
############################################################################################################
decisionRule_warning_message_CUSUM <- function(prodata, data.metrics, method, peptideThresholdRed, peptideThresholdYellow,metricThresholdRed,metricThresholdYellow, L,U, type, listMean,listSD, guidset_selected ) {


  a1 <- number.Of.Out.Of.Range.Metrics(prodata,data.metrics,method = "CUSUM", peptideThresholdRed,peptideThresholdYellow,
                                       L, U, type = 1,listMean,listSD, guidset_selected)
  CUSUMCounterAboveRed1 <- a1[1]
  CUSUMCounterAboveYellow1 <- a1[2]

  a2 <- number.Of.Out.Of.Range.Metrics(prodata,data.metrics,method = "CUSUM", peptideThresholdRed,peptideThresholdYellow,
                                       L, U, type = 2,listMean,listSD, guidset_selected)
  CUSUMCounterAboveRed2 <- a2[1]
  CUSUMCounterAboveYellow2 <-a2[2]

  if(CUSUMCounterAboveRed1 > metricThresholdRed && CUSUMCounterAboveRed2 > metricThresholdRed) return({paste("CUSUM RED FLAG: System performance is UNACCEPTABLE")})
  if(CUSUMCounterAboveRed1 > metricThresholdRed && CUSUMCounterAboveRed2 <= metricThresholdRed) return({paste("CUSUM RED FLAG: System performance is UNACCEPTABLE")})
  if(CUSUMCounterAboveRed1 <= metricThresholdRed && CUSUMCounterAboveRed2 > metricThresholdRed) return({paste("CUSUM RED FLAG: System performance is UNACCEPTABLE")})
  if(CUSUMCounterAboveYellow1 > metricThresholdYellow && CUSUMCounterAboveYellow2 > metricThresholdYellow) return({"CUSUM Yellow FLAG: System performance is POOR"})
  if(CUSUMCounterAboveYellow1 <= metricThresholdYellow && CUSUMCounterAboveYellow2 > metricThresholdYellow) return({"CUSUM Yellow FLAG: System performance is POOR"})
  if(CUSUMCounterAboveYellow1 > metricThresholdYellow && CUSUMCounterAboveYellow2 <= metricThresholdYellow) return({"CUSUM Yellow FLAG: System performance is POOR"})

  return({"CUSUM BLUE FLAG: System performance is acceptable"})

  }
