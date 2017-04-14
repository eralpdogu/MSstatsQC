getMetricData <- function(prodata, peptide, L, U, metric, normalization, selectMean, selectSD) {
  #"Precursor" is one of the columns in data that shows the name of peptides
  precursor.data<-prodata[prodata$Precursor==peptide,]
  metricData <- 0
  mu <- 0
  sd <- 0

  if(is.null(metric)){
    return(NULL)
  }

  metricData = precursor.data[,metric]

  if(normalization == TRUE) {
    if(is.null(selectMean) && is.null(selectSD)) {
      mu=mean(metricData[L:U]) # in-control process mean
      sd=sd(metricData[L:U]) # in-control process variance
      # print(mu)
      # print(sd)
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

  # if(normalization == TRUE) {
  #   mu=mean(metricData[L:U]) # in-control process mean
  #   sd=sd(metricData[L:U]) # in-control process variance
  #   if(sd == 0) {sd <- 0.0001}
  #   metricData=scale(metricData[1:length(metricData)],mu,sd) # transformation for N(0,1) )
  #   return(metricData)
  # } else if(normalization == FALSE){
  #   return(metricData)
  # }

}
