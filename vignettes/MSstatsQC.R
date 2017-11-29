## ----eval=TRUE-----------------------------------------------------------
source("https://bioconductor.org/biocLite.R")
biocLite("MSstatsQC")

## ---- eval=TRUE----------------------------------------------------------
#A typical multi peptide and multi metric system suitability dataset
#This dataset was generated during CPTAC Study 9.1 at Site 54
library(MSstatsQC)
data <- MSstatsQC::S9Site54

## ---- eval=FALSE---------------------------------------------------------
#  MSnbaseToMSstatsQC(msfile)

## ---- eval=TRUE----------------------------------------------------------
data<-DataProcess(data)

## ---- eval=TRUE----------------------------------------------------------
#An X chart when a guide set (1-20 runs) is used to monitor the mean of retention time
XmRChart( data, peptide = "TAAYVNAIEK", L = 1, U = 20, metric = "BestRetentionTime", normalization = FALSE, ytitle = "X Chart : retention time", type = "mean", selectMean = NULL ,selectSD = NULL )
#An X chart when a guide set (1-20 runs) is used to monitor the mean of total peak area
XmRChart( data, peptide = "TAAYVNAIEK", L = 1, U = 20, metric = "TotalArea", normalization = FALSE, ytitle = "X Chart : peak area", type = "mean", selectMean = NULL ,selectSD = NULL )
#An X chart when a guide set (1-20 runs) is used to monitor the variability of retention time
XmRChart( data, peptide = "TAAYVNAIEK", L = 1, U = 20, metric = "BestRetentionTime", normalization = FALSE, ytitle = "mR Chart : retention time", type = "variability", selectMean = NULL ,selectSD = NULL )
#An X chart when a guide set (1-20 runs) is used to monitor the variability of total peak area
XmRChart( data, peptide = "TAAYVNAIEK", L = 1, U = 20, metric = "TotalArea", normalization = FALSE, ytitle = "mR Chart : peak area", type = "variability", selectMean = NULL, selectSD = NULL )
#An X chart when a guide set (1-20 runs) is used to monitor the mean of retention time
XmRChart( data, peptide = "TAAYVNAIEK", L = 1, U = 20, metric = "BestRetentionTime", normalization = FALSE, ytitle = "X Chart : retention time", type = "mean", selectMean = NULL, selectSD = NULL )
#An X chart when a guide set (1-20 runs) is used to monitor the mean of total peak area
XmRChart( data, peptide = "TAAYVNAIEK", L = 1, U = 20, metric = "TotalArea", normalization = FALSE, ytitle = "X Chart : peak area", type = "mean", selectMean = NULL, selectSD = NULL )
#An X chart when a guide set (1-20 runs) is used to monitor the variability of retention time
XmRChart( data, peptide = "TAAYVNAIEK", L = 1, U = 20, metric = "BestRetentionTime", normalization = FALSE, ytitle = "mR Chart : retention time", type = "variability", selectMean = NULL, selectSD = NULL )
#An X chart when a guide set (1-20 runs) is used to monitor the variability of total peak area
XmRChart( data, peptide = "TAAYVNAIEK", L = 1, U = 20, metric = "TotalArea", normalization = FALSE, ytitle = "mR Chart : peak area", type = "variability", selectMean = NULL, selectSD = NULL )
#Mean and standard deviation of LVNELTEFAK is known
XmRChart( data, "LVNELTEFAK", metric = "BestRetentionTime", selectMean = 28.5, selectSD = 1 )

## ---- eval=TRUE, echo =FALSE, fig.height=3-------------------------------
#A CUSUMm chart when a guide set (1-20 runs) is used to monitor the mean of retention time
CUSUMChart( data, peptide = "TAAYVNAIEK", L = 1, U = 20, metric = "BestRetentionTime", normalization = TRUE, ytitle = "CUSUMm Chart : retention time", type = "mean", referenceValue = 0.5, decisionInterval = 5, selectMean = NULL ,selectSD = NULL )
#A CUSUMm chart when a guide set (1-20 runs) is used to monitor the mean of total peak area
CUSUMChart( data, peptide = "TAAYVNAIEK", L = 1, U = 20, metric = "TotalArea", normalization = TRUE, ytitle = "CUSUMm Chart : peak area", type = "mean", referenceValue = 0.5, decisionInterval = 5, selectMean = NULL ,selectSD = NULL  )
#A CUSUMv chart when a guide set (1-20 runs) is used to monitor the variability of retention time
CUSUMChart( data, peptide = "TAAYVNAIEK", L = 1, U = 20, metric = "BestRetentionTime", normalization = TRUE, ytitle = "CUSUMv Chart : retention time", type = "variability", referenceValue = 0.5, decisionInterval = 5, selectMean = NULL ,selectSD = NULL )
#A CUSUMv chart when a guide set (1-20 runs) is used to monitor the variability of total peak area
CUSUMChart( data, peptide = "TAAYVNAIEK", L = 1, U = 20, metric = "TotalArea", normalization = TRUE, ytitle = "CUSUMv Chart : peak area", type = "variability", referenceValue = 0.5, decisionInterval = 5, selectMean = NULL ,selectSD = NULL )
#A CUSUMm chart when a guide set (1-20 runs) is used to monitor the mean of retention time
CUSUMChart( data, peptide = "TAAYVNAIEK", L = 1, U = 20, metric = "BestRetentionTime", normalization = TRUE, ytitle = "CUSUMm Chart : retention time", type = "mean", referenceValue = 0.5, decisionInterval = 5, selectMean = NULL ,selectSD = NULL  )
#A CUSUMm chart when a guide set (1-20 runs) is used to monitor the mean of total peak area
CUSUMChart( data, peptide = "TAAYVNAIEK", L = 1, U = 20, metric = "TotalArea", normalization = TRUE, ytitle = "CUSUMm Chart : peak area", type = "mean", referenceValue = 0.5, decisionInterval = 5, selectMean = NULL ,selectSD = NULL  )
#A CUSUMv chart when a guide set (1-20 runs) is used to monitor the variability of retention time
CUSUMChart( data, peptide = "TAAYVNAIEK", L = 1, U = 20, metric = "BestRetentionTime", normalization = TRUE, ytitle = "mR Chart : retention time", type = "variability", referenceValue = 0.5, decisionInterval = 5, selectMean = NULL ,selectSD = NULL )
#A CUSUMv chart when a guide set (1-20 runs) is used to monitor the variability of total peak area
CUSUMChart( data, peptide = "TAAYVNAIEK", L = 1, U = 20, metric = "TotalArea", normalization = TRUE, ytitle = "CUSUMv Chart : peak area", type = "variability", referenceValue = 0.5, decisionInterval = 5, selectMean = NULL ,selectSD = NULL )

## ---- eval=FALSE---------------------------------------------------------
#  # Retention time >> first 20 observations are used as a guide set
#  XmRChart(data, "TAAYVNAIEK", metric = "BestRetentionTime", type="mean", L = 1, U = 20)
#  ChangePointEstimator(data, "TAAYVNAIEK", metric = "BestRetentionTime", type="mean", L = 1, U = 20)

## ---- eval=TRUE----------------------------------------------------------
# Retention time >> first 20 observations are used as a guide set
XmRChart(data, "YSTDVSVDEVK", metric = "BestRetentionTime", type="mean", L = 1, U = 20)
ChangePointEstimator(data, "YSTDVSVDEVK", metric = "BestRetentionTime", type="variability", L = 1, U = 20)

## ---- eval=TRUE, echo=FALSE----------------------------------------------
# Retention time >> first 20 observations are used as a guide set
RiverPlot(data = S9Site54, L = 1, U = 20, method = "XmR")
RiverPlot(data = S9Site54, L = 1, U = 20, method = "CUSUM")
RadarPlot(data = S9Site54, L = 1, U = 20, method = "XmR")
RadarPlot(data = S9Site54, L = 1, U = 20, method = "CUSUM")

## ---- eval=TRUE, echo=TRUE-----------------------------------------------
# A decision map for Site 54 can be generated using the following script
# Retention time >> first 20 observations are used as a guide set
DecisionMap(data,method="XmR",peptideThresholdRed = 0.25,peptideThresholdYellow = 0.10,
                         L = 1, U = 20,type = "mean",title = "Decision map",listMean = NULL,listSD = NULL)

## ---- eval=FALSE---------------------------------------------------------
#  #Saving plots generated by plotly
#  p<-XmRChart( data, peptide = "TAAYVNAIEK", L = 1, U = 20, metric = "BestRetentionTime", normalization = FALSE,
#                        ytitle = "X Chart : retention time", type = "mean", selectMean = NULL, selectSD = NULL )
#  htmlwidgets::saveWidget(p, "Aplot.html")
#  export(p, file = "Aplot.png")
#  
#  #Saving plots generated by ggplot2
#  p<-RiverPlot(data, L=1, U=20)
#  ggsave(filename="Summary.pdf", plot=p)
#  #or
#  ggsave(filename="Summary.png", plot=p)

