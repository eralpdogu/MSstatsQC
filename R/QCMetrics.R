#INPUTS : "prodata" is the data user uploads.
#         "precursorSelection" is the precursor that user selects in Data Import tab.
#                              it can be either one precursor(peptide) or it can be "all peptides"
#         "L" and "U" are lower and upper bound of guide set that user choose in Data Import tab.
#         "metric" is a metric that user defines in his data set
#         "normalization" is either TRUE or FALSE
#Description of function : it gets the metric column only for the precursor chosen
#and either return the column as it is or normalize it and then return it
# getMetricData <- function(prodata, precursorSelection, L, U, metric, normalization) {
#   precursor.data<-prodata[prodata$Precursor==precursorSelection,]
#   metricData <- 0
#
#   if(is.null(metric)){
#     return(NULL)
#   }
#
#   metricData = precursor.data[,metric]
#   if(normalization == TRUE) {
#     mu=mean(metricData[L:U]) # in-control process mean
#     sd=sd(metricData[L:U]) # in-control process variance
#     if(sd == 0) {sd <- 0.0001}
#     metricData=scale(metricData[1:length(metricData)],mu,sd) # transformation for N(0,1) )
#     return(metricData)
#   } else if(normalization == FALSE){
#     return(metricData)
#   }
#
# }
#########################################################################################################
find_custom_metrics <- function(prodata) {

  prodata <- prodata[, which(colnames(prodata)=="Annotations"):ncol(prodata),drop = FALSE]
  nums <- sapply(prodata, is.numeric)
  # limiting custom metrics up to 10 metrics and not more
  other.metrics <- colnames(prodata[,nums])[1:ifelse(length(colnames(prodata[,nums]))<11,
                                                     length(colnames(prodata[,nums])),
                                                     10)
                                            ]
  if(any(is.na(other.metrics))) {
    return(c())
  }
  return(other.metrics)
}
################################################################
#z = metricData
#INPUT : "prodata" is the data user uploads.
#        "metricData" is the column of the data related to the metric we want.
#                     Forexample if we want retention time, it gives retention time column
#        "precursorSelection" is the precursor that user selects in Data Import tab.
#                             it can be either one precursor(peptide) or it can be "all peptides"
#        "L" and "U" are lower and upper bound of guide set that user choose in Data Import tab.
#        "type" is either "mean" or "dispersion". one is "Individual Value" plot and other
#               "Moving Range" plot
#DESCRIPTION : returns a data frame for CUSUM that contains all the information needed to plot CUSUM
CUSUM.data.prepare <- function(prodata, metricData, precursorSelection, type) {

  k=0.5
  CUSUM.outrange.thld <- 5
  outRangeInRangePoz = rep(0,length(metricData))
  outRangeInRangeNeg = rep(0,length(metricData))
  #precursor.data only gets the data for the selected precursor
  #"Precursor" is one of the columns in data that shows the name of peptides
  precursor.data <- prodata[prodata$Precursor==precursorSelection,]

  v <- numeric(length(metricData))

  Cpoz <- numeric(length(metricData))
  Cneg <- numeric(length(metricData))

  for(i in 2:length(metricData)) {
    Cpoz[i]=max(0,(metricData[i]-(k)+Cpoz[i-1]))
    Cneg[i]=max(0,((-k)-metricData[i]+Cneg[i-1]))
  }

  if(type == "dispersion") {
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
#        "metricData" is the column of the data related to the metric we want.
#         Forexample if we want retention time, it gives retention time column
#        "type" is either "mean" or "dispersion". one is "Individual Value" plot and
#               other "Moving Range" plot
#DESCRIPTION : returns a data frame for CP that contains all the information needed to plot
#Change Point
CP.data.prepare <- function(prodata, metricData, type) {

  Et <-  numeric(length(metricData)-1) # this is Ct in type = mean , and Dt in type = dispersion.
  SS<- numeric(length(metricData)-1)
  SST<- numeric(length(metricData)-1)
  tho.hat <- 0

  if(type == "mean") {
    ## Change point analysis for mean (Single step change model)
    for(i in 1:(length(metricData)-1)) {
      Et[i]=(length(metricData)-i)*
        (((1/(length(metricData)-i))*
            sum(metricData[(i+1):length(metricData)]))-0)^2 #change point function
    }
    QCno=1:(length(metricData)-1)
  } else if(type == "dispersion") {
    ## Change point analysis for variance (Single step change model)
    for(i in 1:length(metricData)) {
      SS[i]=metricData[i]^2
    }
    for(i in 1:length(metricData)) {
      SST[i]=sum(SS[i:length(metricData)])
      #change point function
      Et[i]=((SST[i]/2)-((length(metricData)-i+1)/2)*
               log(SST[i]/(length(metricData)-i+1))-(length(metricData)-i+1)/2)
    }
    QCno=1:length(metricData)
  }
  tho.hat = which(Et==max(Et)) # change point estimate
  plot.data <- data.frame(QCno,Et,tho.hat)

  return(plot.data)
}
###################################################################################################
#INPUT : "prodata" is the data user uploads.
#        "L" and "U" are lower and upper bound of guide set that user choose in Data Import tab.
#        "data.metrics" is all the available metrics. It is defined in server.R
get_CP_tho.hat <- function(prodata, L, U, data.metrics, listMean, listSD) {
  tho.hat <- data.frame(tho.hat = c(), metric = c(), group = c(), y=c())
  precursors <- levels(prodata$Precursor)
  for(metric in data.metrics) {
    for (j in 1:nlevels(prodata$Precursor)) {
      metricData <- getMetricData(prodata, precursors[j], L, U, metric = metric,normalization = TRUE,
                                  selectMean = listMean[[metric]],selectSD = listSD[[metric]])
      mix <- rbind(
        data.frame(tho.hat = CP.data.prepare(prodata, metricData, type = "mean")$tho.hat[1],
                   metric = metric, group = "Individual Value", y=1.1),
        data.frame(tho.hat = CP.data.prepare(prodata, metricData, type = "dispersion")$tho.hat[1],
                   metric = metric, group = "Moving Range", y=-1.1)
      )
      tho.hat <- rbind(tho.hat, mix)

    }
  }

  return(tho.hat)
}
###################################################################################################
#z = metricData
#INPUT : "prodata" is the data user uploads.
#        "metricData" is the column of the data related to the metric we want.
#                     Forexample if we want retention time, it gives retention time column
#        "L" and "U" are lower and upper bound of guide set that user choose in Data Import tab.
#        "type" is either "mean" or "dispersion". one is "Individual Value" plot and
#         other "Moving Range" plot
#DESCRIPTION : returns a data frame for XmR that contains all the information needed to plot XmR
XmR.data.prepare <- function(prodata, metricData, L,U, type, selectMean,selectSD) {
  t <- numeric(length(metricData)-1)
  UCL <- 0
  LCL <- 0
  InRangeOutRange <- rep(0,length(metricData))

  for(i in 2:length(metricData)) {
    t[i]=abs(metricData[i]-metricData[i-1]) # Compute moving range of metricData
  }

  QCno=1:length(metricData)

  if(type == "mean") {
    if(is.null(selectMean) && is.null(selectSD)) {
      UCL=mean(metricData[L:U])+2.66*sd(t[L:U])
      LCL=mean(metricData[L:U])-2.66*sd(t[L:U])
    }else {
      UCL = selectMean + 2.66 * selectSD
      LCL = selectMean - 2.66 * selectSD
    }
    t <- metricData
  }else if(type == "dispersion") {
    ## Calculate MR chart statistics and limits
    if(is.null(selectMean) && is.null(selectSD)) {
      UCL=3.267*sd(t[1:L-U])
    }else{
      UCL = 3.267 * selectSD
    }
    LCL=0
  }
  # if(type == "mean") {
  #     UCL=mean(metricData[L:U])+2.66*sd(t[L:U])
  #     LCL=mean(metricData[L:U])-2.66*sd(t[L:U])
  #     t <- metricData
  #
  # } else if(type == "dispersion") {
  ## Calculate MR chart statistics and limits
  #     UCL=3.267*sd(t[1:L-U])
  #     LCL=0
  # }

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
#         "metric" is a metric that user defines in his data set
#         "L" and "U" are lower and upper bound of guide set that user choose in Data Import tab.
#        "type" is either "mean" or "dispersion". one is "Individual Value" plot and
#               other "Moving Range" plot
#DESCRIPTION : returns a data frame that is used in CUSUM.Summary.DataFrame function below to
#              use for summary plot of CUSUM in Summary Tab of shiny app
CUSUM.Summary.prepare <- function(prodata, metric, L, U,type, selectMean, selectSD) {
  h <- 5

  QCno <- 1:nrow(prodata)
  y.poz <- rep(0,nrow(prodata))
  y.neg <- rep(0,nrow(prodata))
  counter <- rep(0,nrow(prodata))

  precursors <- levels(prodata$Precursor)

  for(j in 1:length(precursors)) {
    metricData <- getMetricData(prodata, precursors[j], L, U, metric = metric,
                                normalization = TRUE, selectMean, selectSD)
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
                          group = ifelse(rep(type == "mean",2*max_QCno),
                                         c(rep("Mean increase",max_QCno),
                                           rep("Mean decrease",max_QCno)),
                                         c(rep("Variability increase",max_QCno),
                                           rep("Variability decrease",max_QCno))),
                          metric = rep(metric,max_QCno*2)
  )
  return(plot.data)
}
############################################################################################
CUSUM.Summary.DataFrame <- function(prodata, data.metrics, L, U,listMean,listSD) {
  dat <- data.frame(QCno = c(),
                    pr.y = c(),
                    group = c(),
                    metric = c())
  for (metric in data.metrics) {
    data.1   <- CUSUM.Summary.prepare(prodata, metric = metric, L, U,type = "mean",
                                      selectMean = listMean[[metric]], selectSD = listSD[[metric]])
    data.2   <- CUSUM.Summary.prepare(prodata, metric = metric, L, U,type = "dispersion",
                                      selectMean = listMean[[metric]], selectSD = listSD[[metric]])
    data.2$pr.y <- -(data.2$pr.y)
    dat <- rbind(dat,data.1,data.2)
  }
  return(dat)
}
############################################################################################
#DESCRIPTION : for each metric returns a data frame of QCno, probability of out of control
#peptide for dispersion or mean plot
XmR.Summary.prepare <- function(prodata, metric, L, U,type,selectMean,selectSD) {
  QCno    <- 1:nrow(prodata)
  y.poz <- rep(0,nrow(prodata))
  y.neg <- rep(0,nrow(prodata))
  counter <- rep(0,nrow(prodata))

  precursors <- levels(prodata$Precursor)
  #print(selectMean)
  for(j in 1:length(precursors)) {
    metricData <- getMetricData(prodata, precursors[j], L = L, U = U, metric = metric,
                                normalization = TRUE,selectMean,selectSD)
    counter[1:length(metricData)] <- counter[1:length(metricData)]+1
    plot.data <- XmR.data.prepare(prodata, metricData , L , U , type,selectMean,selectSD)

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
                          group = ifelse(rep(type == "mean",2*max_QCno),
                                         c(rep("Mean increase",max_QCno),
                                           rep("Mean decrease",max_QCno)),
                                         c(rep("Variability increase",max_QCno),
                                           rep("Variability decrease",max_QCno))),
                          metric = rep(metric,max_QCno*2))

  return(plot.data)
}
###########################################################################################
#DESCRIPTION : returns a data frame for all the metrics, probability of out of range peptides,
#wheter it is for metric mean or metric dispersion and i is an increase or decrease
XmR.Summary.DataFrame <- function(prodata, data.metrics, L, U, listMean, listSD) {
  dat <- data.frame(QCno = c(),
                    pr.y = c(),
                    group = c(),
                    metric = c())

  for (metric in data.metrics) {
    #print("in DF")
    #print(listMean[[metric]])
    data.1   <- XmR.Summary.prepare(prodata, metric = metric, L, U,type = "mean",
                                    selectMean = listMean[[metric]],selectSD = listSD[[metric]])
    data.2   <- XmR.Summary.prepare(prodata, metric = metric, L, U,type = "dispersion",
                                    selectMean = listMean[[metric]],selectSD = listSD[[metric]])
    data.2$pr.y <- -(data.2$pr.y)
    dat <- rbind(dat, data.1, data.2)
  }
  return(dat)
}
############################################################################################
heatmap.DataFrame <- function(prodata, data.metrics,method,
                              peptideThresholdRed,peptideThresholdYellow, L, U,
                              type, listMean, listSD) {
  #peptideThresholdGood = peptideThresholdRed
  #peptideThresholdWarn = peptideThresholdYellow
  time <- c()
  val <- c()
  met <- c()
  bin <- c()

  for (metric in data.metrics) {
    df <- Decision.DataFrame.prepare(prodata, metric, method,
                                     peptideThresholdRed,peptideThresholdYellow,
                                     L, U,type, selectMean = listMean[[metric]],
                                     selectSD=listSD[[metric]])
    #print(df)
    time_df <- as.character(df$AcquiredTime)
    val_df <- df$pr.y
    met_df <- rep(metric,length(val_df))
    bin_df <- df$bin
    time <- c(time,time_df)
    val <- c(val,val_df)
    met <- c(met,met_df)
    bin <- c(bin,bin_df)
  }

  dataFrame <- data.frame(time = time,
                          value = val,
                          metric = met,
                          bin = bin
  )
  return(dataFrame)
}
############################################################################################
Compute.QCno.OutOfRangePeptide.XmR <- function(prodata,L,U,metric,type,
                                               XmR.type,selectMean,selectSD) {
  precursors <- levels(prodata$Precursor)
  QCno.out.range <- c()

  for(j in 1:length(precursors)) {
    metricData <- getMetricData(prodata, precursors[j], L = L, U = U,
                                metric = metric, normalization = TRUE,selectMean,selectSD)
    plot.data <- XmR.data.prepare(prodata, metricData , L = L, U = U, type,selectMean,selectSD)
    if(XmR.type == "poz")
      QCno.out.range <- c(QCno.out.range,length(plot.data[plot.data$t >= plot.data$UCL, ]$QCno))
    else
      QCno.out.range <- c(QCno.out.range,length(plot.data[plot.data$t <= plot.data$LCL, ]$QCno))
  }
  return(QCno.out.range)
}
#############################################################################################
Compute.QCno.OutOfRangePeptide.CUSUM <- function(prodata,L,U,metric,type,
                                                 CUSUM.type,selectMean,selectSD) {
  h <- 5
  precursors <- levels(prodata$Precursor)
  QCno.out.range <- c()

  for(j in 1:length(precursors)) {
    metricData <- getMetricData(prodata, precursors[j], L, U, metric = metric,
                                normalization = TRUE,selectMean,selectSD)
    plot.data <- CUSUM.data.prepare(prodata, metricData, precursors[j], type)
    if(CUSUM.type == "poz")
      QCno.out.range <- c(QCno.out.range,
                          length(plot.data[plot.data$CUSUM.poz >= h | plot.data$CUSUM.poz <= -h, ]$QCno))
    else
      QCno.out.range <- c(QCno.out.range,
                          length(plot.data[plot.data$CUSUM.neg >= h | plot.data$CUSUM.neg <= -h, ]$QCno))
  }
  return(QCno.out.range)
}
###############################################################################################################
XmR.Radar.Plot.prepare <- function(prodata,L,U, metric, type,group,
                                   XmR.type,selectMean,selectSD) {
  precursors <- levels(prodata$Precursor)
  precursors2 <- substring(precursors, first = 1, last = 3)
  QCno.length <- c()
  QCno.out.range.poz <- c() #add
  QCno.out.range.neg <- c() #add
  for(j in 1:length(precursors)) {
    metricData <- getMetricData(prodata, precursors[j], L = L, U = U,
                                metric = metric, normalization = TRUE,
                                selectMean,selectSD)
    QCno.length <- c(QCno.length,length(metricData))
    plot.data <- XmR.data.prepare(prodata, metricData , L = L, U = U,
                                  type ,selectMean,selectSD) #add
    QCno.out.range.poz <- c(QCno.out.range.poz,
                            length(plot.data[plot.data$t >= plot.data$UCL, ]$QCno)) #add
    QCno.out.range.neg <- c(QCno.out.range.neg,
                            length(plot.data[plot.data$t <= plot.data$LCL, ]$QCno)) #add
  }
  #add
  if(XmR.type == "poz") {
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
}
################################################################################################
XmR.Radar.Plot.DataFrame <- function(prodata, data.metrics, L,U, listMean, listSD) {
  dat <- data.frame(peptides = c(), OutRangeQCno = c(), group = c(),
                    orderby = c(), metric = c(), tool = c(),
                    probability   = c()
  )
  for (metric in data.metrics) {
    data.1 <- XmR.Radar.Plot.prepare(prodata,L,U,metric = metric,
                                     type = "mean",group = "Mean increase",
                                     XmR.type = "poz",
                                     selectMean = listMean[[metric]],selectSD = listSD[[metric]])
    data.2 <- XmR.Radar.Plot.prepare(prodata,L,U,metric = metric,
                                     type = "mean",group = "Mean decrease",
                                     XmR.type = "neg",
                                     selectMean = listMean[[metric]],selectSD = listSD[[metric]])
    data.3 <- XmR.Radar.Plot.prepare(prodata,L,U,metric = metric,
                                     type = "dispersion",group = "Variability increase",
                                     XmR.type = "poz",
                                     selectMean = listMean[[metric]],selectSD = listSD[[metric]])
    data.4 <- XmR.Radar.Plot.prepare(prodata,L,U,metric = metric,
                                     type = "dispersion",group = "Variability decrease",
                                     XmR.type = "neg",
                                     selectMean = listMean[[metric]],selectSD = listSD[[metric]])
    dat <- rbind(dat, data.1, data.2, data.3, data.4)
  }
  return(dat)
}
#################################################################################################################
CUSUM.Radar.Plot.prepare <- function(prodata,L,U, metric,type,group,
                                     CUSUM.type,selectMean, selectSD) {
  h <- 5
  precursors <- levels(prodata$Precursor)
  precursors2 <- substring(precursors, first = 1, last = 3)
  QCno.length <- c()
  QCno.out.range <- c() #add
  for(j in 1:length(precursors)) {
    metricData <- getMetricData(prodata, precursors[j], L = L, U = U,
                                metric = metric, normalization = TRUE, selectMean, selectSD)
    QCno.length <- c(QCno.length,length(metricData))
    plot.data <- CUSUM.data.prepare(prodata, metricData, precursors[j], type) #add
    #add
    if(CUSUM.type == "poz")
      QCno.out.range <- c(QCno.out.range,
                          length(plot.data[plot.data$CUSUM.poz >= h |
                                           plot.data$CUSUM.poz <= -h, ]$QCno))
    else
      QCno.out.range <- c(QCno.out.range,
                          length(plot.data[plot.data$CUSUM.neg >= h |
                                           plot.data$CUSUM.neg <= -h, ]$QCno))
  }
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
CUSUM.Radar.Plot.DataFrame <- function(prodata, data.metrics, L,U, listMean, listSD) {
  dat <- data.frame(peptides = c(), OutRangeQCno = c(), group = c(),
                    orderby = c(), metric = c(), tool = c(),
                    probability   = c()
  )
  for (metric in data.metrics) {
    data.1 <- CUSUM.Radar.Plot.prepare(prodata,L,U, metric = metric, type = "mean",
                                       group = "Mean increase",CUSUM.type = "poz",
                                       selectMean = listMean[[metric]],selectSD = listSD[[metric]])
    data.2 <- CUSUM.Radar.Plot.prepare(prodata,L,U, metric = metric, type = "mean",
                                       group = "Mean decrease", CUSUM.type = "neg",
                                       selectMean = listMean[[metric]],selectSD = listSD[[metric]])
    data.3 <- CUSUM.Radar.Plot.prepare(prodata,L,U, metric = metric, type = "dispersion",
                                       group = "Variability increase", CUSUM.type = "poz",
                                       selectMean = listMean[[metric]],selectSD = listSD[[metric]])
    data.4 <- CUSUM.Radar.Plot.prepare(prodata,L,U, metric = metric, type = "dispersion",
                                       group = "Variability decrease", CUSUM.type = "neg",
                                       selectMean = listMean[[metric]],selectSD = listSD[[metric]])
    dat <- rbind(dat, data.1, data.2, data.3, data.4)
  }
  return(dat)
}
#######################################################################################################
Decision.DataFrame.prepare <- function(prodata, metric, method, peptideThresholdRed,
                                       peptideThresholdYellow, L, U,type, selectMean, selectSD) {

  h <- 5
  AcquiredTime <- prodata$AcquiredTime
  QCno    <- 1:nrow(prodata)
  y <- rep(0,nrow(prodata))
  counter <- rep(0,nrow(prodata))

  precursors <- levels(prodata$Precursor)
  #add
  if(method == "XmR") {
    for(precursor in precursors) {
      metricData <- getMetricData(prodata, precursor, L = L, U = U,
                                  metric = metric, normalization = TRUE,selectMean,selectSD)
      counter[1:length(metricData)] <- counter[1:length(metricData)]+1
      plot.data <- XmR.data.prepare(prodata, metricData , L , U , type,selectMean,selectSD)
      sub <- plot.data[plot.data$InRangeOutRange == "OutRange",]
      #sub1 <- plot.data[(plot.data$t >= plot.data$UCL), ]
      #sub2 <- plot.data[plot.data$t <= plot.data$LCL, ]
      #sub <- rbind(sub1,sub2)
      y[sub$QCno] <- y[sub$QCno] + 1
    }
  } else if(method == "CUSUM") {
    for(precursor in precursors) {
      metricData <- getMetricData(prodata, precursor, L = L, U = U,
                                  metric = metric, normalization = TRUE,
                                  selectMean,selectSD)
      counter[1:length(metricData)] <- counter[1:length(metricData)]+1
      plot.data <- CUSUM.data.prepare(prodata, metricData, precursor, type)
      sub <- plot.data[(plot.data$CUSUM.poz >= h |
                          plot.data$CUSUM.poz <= -h) |
                         (plot.data$CUSUM.neg >= h |
                            plot.data$CUSUM.neg <= -h), ]
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
                          bin = rep(0,max_QCno)
  )



  for (i in 1:max_QCno) {
    if(plot.data$pr.y[i] > peptideThresholdRed){
      plot.data$bin[i] <- "Fail"
    }
    else if(plot.data$pr.y[i] > peptideThresholdYellow){
      plot.data$bin[i] <- "Warning"
    }
    else {
      plot.data$bin[i] <- "Pass"
    }
  }

  if(type == 2) {

    return(plot.data[-1,])
  }
  return(plot.data)
}
#######################################################################################################
number.Of.Out.Of.Range.Metrics <- function(prodata,data.metrics,method,
                                           peptideThresholdRed, peptideThresholdYellow,
                                           L, U, type, listMean, listSD) {

  metricCounterAboveRed = 0
  metricCounterAboveYellowBelowRed = 0
  precursors <- levels(prodata$Precursor) #add
  #add
  for (metric in data.metrics) {
    QCno    <- 1:nrow(prodata)
    y <- rep(0,nrow(prodata))
    counter <- rep(0,nrow(prodata))
    for(precursor in precursors) {
      metricData <- getMetricData(prodata, precursor, L = L, U = U, metric = metric,
                                  normalization = TRUE,listMean[[metric]],listSD[[metric]])
      counter[1:length(metricData)] <- counter[1:length(metricData)]+1
      #if(method == "XmR") {
      plot.data <- XmR.data.prepare(prodata, metricData , L , U ,type,
                                    selectMean = listMean[[metric]],selectSD = listSD[[metric]])
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
