require(testthat)
require(MSstatsQC)
###########################################################################################
test_that("X and mR functions are calculated correctly",{
  metricData <- getMetricData(data = S9Site54, peptide = "VLVLDTDYK", 
                              L = 1, U = 10,metric = "BestRetentionTime", 
                              normalization = FALSE, selectMean = NULL, selectSD = NULL)
  
  df_func <- XmR.data.prepare(metricData = metricData, L=1, U=10, type = "mean", selectMean = NULL, selectSD = NULL)[1:5,]
  
  df_func$IndividualValue <- floor(df_func$IndividualValue)
  df_func$mR <- floor(df_func$mR)
  df_func$UCL <- floor(df_func$UCL)
  df_func$LCL <- floor(df_func$LCL)

  a <-  c(1,1,1,1,1)
  InRangeOutRange <- factor(a,levels=1:2)
  levels(InRangeOutRange) <- c("InRange","OutRange")
  df <- data.frame(QCno = seq(1,5), IndividualValue = c(24,24,24,24,24),
                   mR = c(24,24,24,24,24), UCL=c(24,24,24,24,24), LCL=c(24,24,24,24,24), 
                   InRangeOutRange)
  
  expect_equal(df_func,df)
})
#############################################################################################
test_that("CUSUM functions are calculated correctly for changes in variability of a metric",{
  metricData <- getMetricData(data = S9Site54, peptide = "VLVLDTDYK", 
                              L = 1, U = 5,metric = "BestRetentionTime", 
                              normalization = FALSE, selectMean = NULL, selectSD = NULL)
  
  df_func <- CUSUM.data.prepare(data = S9Site54, metricData = metricData, 
                                peptide = "VLVLDTDYK", type = "variability", 
                                referenceValue = 0.5, decisionInterval = 5)[1:5,]
  
  df_func$Annotations <- NA
  df_func$CUSUM.poz <- floor(df_func$CUSUM.poz)
  df_func$CUSUM.neg <- floor(df_func$CUSUM.neg)
  
  df <- data.frame(QCno = seq(1,5), CUSUM.poz = c(0,11,22,34,45),
                   CUSUM.neg = c(0,0,0,0,0), Annotations = c(1,1,1,1,1),
                   outRangeInRangePoz = c("InRangeCUSUM+","OutRangeCUSUM+","OutRangeCUSUM+","OutRangeCUSUM+","OutRangeCUSUM+"),
                   outRangeInRangeNeg = c("InRangeCUSUM-","InRangeCUSUM-","InRangeCUSUM-","InRangeCUSUM-","InRangeCUSUM-"))
  
  df$Annotations <- NA
  
  expect_equal(df_func,df)
})
###########################################################################################
test_that("CUSUM functions are calculated correctly for changes in mean of a metric",{
metricData <- getMetricData(data = S9Site54, peptide = "VLVLDTDYK", 
                              L = 1, U = 5,metric = "BestRetentionTime", 
                              normalization = FALSE, selectMean = NULL, selectSD = NULL)
  
df_func <- CUSUM.data.prepare(data = S9Site54, metricData = metricData, 
                                peptide = "VLVLDTDYK", type = "mean", 
                                referenceValue = 0.5, decisionInterval = 5)[1:5,]
  
df_func$Annotations <- NA
df_func$CUSUM.poz <- floor(df_func$CUSUM.poz)
df_func$CUSUM.neg <- floor(df_func$CUSUM.neg)

df <- data.frame(QCno = seq(1,5), CUSUM.poz = c(0,24,48,72,96),
                   CUSUM.neg = c(0,0,0,0,0), Annotations = c(1,1,1,1,1),
                   outRangeInRangePoz = c("InRangeCUSUM+","OutRangeCUSUM+","OutRangeCUSUM+","OutRangeCUSUM+","OutRangeCUSUM+"),
                   outRangeInRangeNeg = c("InRangeCUSUM-","InRangeCUSUM-","InRangeCUSUM-","InRangeCUSUM-","InRangeCUSUM-"))
  
df$Annotations <- NA
  
expect_equal(df_func,df)
})
#############################################################################################
test_that("CUSUM functions are calculated correctly for changes in variability of a metric",{
  metricData <- getMetricData(data = S9Site54, peptide = "VLVLDTDYK", 
                              L = 1, U = 5,metric = "BestRetentionTime", 
                              normalization = FALSE, selectMean = NULL, selectSD = NULL)
  
  df_func <- CUSUM.data.prepare(data = S9Site54, metricData = metricData, 
                                peptide = "VLVLDTDYK", type = "variability", 
                                referenceValue = 0.5, decisionInterval = 5)[1:5,]
  
  df_func$Annotations <- NA
  df_func$CUSUM.poz <- floor(df_func$CUSUM.poz)
  df_func$CUSUM.neg <- floor(df_func$CUSUM.neg)
  
  df <- data.frame(QCno = seq(1,5), CUSUM.poz = c(0,11,22,34,45),
                   CUSUM.neg = c(0,0,0,0,0), Annotations = c(1,1,1,1,1),
                   outRangeInRangePoz = c("InRangeCUSUM+","OutRangeCUSUM+","OutRangeCUSUM+","OutRangeCUSUM+","OutRangeCUSUM+"),
                   outRangeInRangeNeg = c("InRangeCUSUM-","InRangeCUSUM-","InRangeCUSUM-","InRangeCUSUM-","InRangeCUSUM-"))
  
  df$Annotations <- NA
  
  expect_equal(df_func,df)
})
#############################################################################################
test_that("Change point functions are calculated correctly for changes in mean of a metric",{
  
metricData <- getMetricData(data = S9Site54, peptide = "VLVLDTDYK", L = 1, U = 5,
                              metric = "BestRetentionTime", normalization = FALSE,
                              selectMean = NULL, selectSD = NULL)
  
df_func <- CP.data.prepare( metricData = metricData, type = "mean")[1:5,]
  
df_func$Et <- floor(df_func$Et)

df <- data.frame(QCno = seq(1,5), Et = c(27444,26834,26232,25627,25021),
                   tho.hat = c(1,1,1,1,1))
  
expect_equal(df_func,df)
})
##############################################################################################
test_that("Change point functions are calculated correctly for changes in variability of a metric",{
  
  metricData <- getMetricData(data = S9Site54, peptide = "VLVLDTDYK", L = 1, U = 5,
                              metric = "BestRetentionTime", normalization = FALSE,
                              selectMean = NULL, selectSD = NULL)
  
  df_func <- CP.data.prepare( metricData = metricData, type = "variability")[1:5,]
  
  df_func$Et <- floor(df_func$Et)
  
  df <- data.frame(QCno = seq(1,5), Et = c(13854,13858,13861,13864,13867),
                   tho.hat = c(46,46,46,46,46))
  
  expect_equal(df_func,df)
})
##############################################################################################
test_that("get_CP_tho.hat first ten rows for reten time works well",{
  df <- c(get_CP_tho.hat(data = S9Site54, L = 1, U = 5, data.metrics = c("BestRetentionTime"), NULL,NULL)[1:10,1])
  expect_equal(df, c(5,45,5,45,5,45,38,45,16,45))
})
