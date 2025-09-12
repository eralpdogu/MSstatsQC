#' A function to train random forest classifiers for QC data
#'
#' @param guide.set comma-separated (.csv), metric file. It should contain a "Precursor" column and the metrics columns. It should also include "Annotations" for each run.
#' @param sim.size enter simulation size.
#' @param guide.set.annotations comma-separated (.csv), metric file with annotations such as pass and fail.
#' @param nfolds fold for cross validation
#' @param a lower threshold to define shift size
#' @param b upper threshold to define shift size
#' @export
#' @import dplyr
#' @import h2o
#' @import ggplot2
#' @import gridExtra
#' @import FrF2
#' @import car
#' @import reshape2
#' @import dplyr
#' @import plyr
#' @return A trained model and performance indicators from train/validation/test splits
#' @keywords machine learning, training, feature engineering
#' @examples
#' S9Site54.dataML<-DataProcess(MSstatsQC::S9Site54[,])
#' colnames(S9Site54.dataML)[1]<-c("idfile")
#' colnames(S9Site54.dataML)[2]<-c("peptide")
#' S9Site54.dataML$peptide <-as.factor(S9Site54.dataML$peptide)
#' S9Site54.dataML$idfile<-as.numeric(S9Site54.dataML$idfile)
#' S9Site54.dataML <- within(S9Site54.dataML, rm(Annotations, missing))
#' guide.set <- filter(S9Site54.dataML, idfile<=20)
#' MSstatsQC.ML.trainR(guide.set, sim.size=10)
#'
MSstatsQC.ML.trainR<- function(guide.set,
                               sim.size,
                               guide.set.annotations=NULL,
                               nfolds=NULL, a = 1.5, b = 2){
  nmetric<-ncol(guide.set)-2
  factor.names = colnames(guide.set[,3:ncol(guide.set)])

  if (is.null(dim(guide.set.annotations))==TRUE) {
    d<-rbind(
    QcClassifier_data_step(guide.set,nmetric,factor.names,sim.size*1,peptide.colname, L=a, U=b),
    QcClassifier_data_var(guide.set,nmetric,factor.names,sim.size*1,peptide.colname, L=a, U=b),
    QcClassifier_data_linear(guide.set,nmetric,factor.names,sim.size*1,peptide.colname, L=a, U=b),

    QcClassifier_data_step(guide.set,nmetric,factor.names,sim.size*1,peptide.colname, L=-b, U=-a),
    QcClassifier_data_var(guide.set,nmetric,factor.names,sim.size*1,peptide.colname, L=-b, U=-a),
    QcClassifier_data_linear(guide.set,nmetric,factor.names,sim.size*1,peptide.colname, L=-b, U=-a)
  )}
  else{
    d<-rbind(
      QcClassifier_data_step(guide.set,nmetric,factor.names,sim.size*1,peptide.colname, L=a, U=b),
      QcClassifier_data_var(guide.set,nmetric,factor.names,sim.size*1,peptide.colname, L=a, U=b),
      QcClassifier_data_linear(guide.set,nmetric,factor.names,sim.size*1,peptide.colname, L=a, U=b),

      QcClassifier_data_step(guide.set,nmetric,factor.names,sim.size*1,peptide.colname, L=-b, U=-a),
      QcClassifier_data_var(guide.set,nmetric,factor.names,sim.size*1,peptide.colname, L=-b, U=-a),
      QcClassifier_data_linear(guide.set,nmetric,factor.names,sim.size*1,peptide.colname, L=-b, U=-a),
    QcClassifier_data_annotated(guide.set, guide.set.annotations)
  )}

  #d<-rbind(d1,d2, d3, d4,d5, d6)
  ## 80% of the sample size
  #smp_size <- floor(0.8 * nrow(d))
  h2o.init()

  d <- as.h2o(d)
  d$RESPONSE <- as.factor(d$RESPONSE)

  #train_ind <- sample(seq_len(nrow(d)), size = smp_size)
  splits <- h2o.splitFrame(
    d,           ##  splitting the H2O frame we read above
    c(0.6,0.2),   ##  create splits of 60% and 20%;
    ##  H2O will create one more split of 1-(sum of these parameters)
    ##  so we will get 0.6 / 0.2 / 1 - (0.6+0.2) = 0.6/0.2/0.2
    seed=1234)    ##  setting a seed will ensure reproducible results (not R's seed)

  train <- h2o.assign(splits[[1]], "train.hex")
  ## assign the first result the R variable train
  ## and the H2O name train.hex
  valid <- h2o.assign(splits[[2]], "valid.hex")   ## R valid, H2O valid.hex
  test <- h2o.assign(splits[[3]], "test.hex")     ## R test, H2O test.hex

  #launch h2o cluster
  localH2O <- h2o.init(nthreads = -1)

  #import r objects to h2o cloud
  train_h2o <- as.h2o(train)
  valid_h2o <- as.h2o(valid)
  test_h2o <- as.h2o(test)

  hyper_grid.h2o <- list(ntrees = seq(50, 200, by = 20),
                         max_depth   = seq(2, 4, by = 1),
                         mtries = seq(3, sqrt(dim(d)[2]), by = 1))

  search_criteria<- list(strategy = "RandomDiscrete",
                            stopping_metric = "AUC",
                            stopping_tolerance = 0.002,
                            stopping_rounds = 2,
                            max_runtime_secs = 30*60)


  if (is.null(nfolds)==TRUE) {
    rf_grid <- h2o.grid(algorithm = "randomForest",
                                              grid_id = "rf_grid",
                                              training_frame = train,        ## the H2O frame for training
                                              validation_frame = valid,      ## the H2O frame for validation (not required)
                                              x= colnames(train_h2o),
                                              y= "RESPONSE",
                                              hyper_params = hyper_grid.h2o)}
  else{
    rf_grid <- h2o.grid(algorithm = "randomForest",
                            grid_id = "rf_grid",
                            training_frame = train,        ## the H2O frame for training
                            validation_frame = valid,      ## the H2O frame for validation (not required)
                            x= colnames(train_h2o),
                            y= "RESPONSE",
                            nfolds=nfolds,
                            hyper_params = hyper_grid.h2o)}
  # Turn parameters for RF:

  grid_perf2 <- h2o.getGrid(grid_id = "rf_grid",
                            sort_by = "AUC",
                            decreasing = FALSE)

  # Best RF:
  rf_model <- h2o.getModel(grid_perf2@model_ids[[1]])

  message(paste("Constructed the RF model"))
  #print(rf_model)
  hit1<-h2o.performance(rf_model, newdata = train)
  hit2<-h2o.performance(rf_model, newdata = valid)
  hit3<-h2o.performance(rf_model, newdata = test)
  message(paste("Performance"))
  print(c(training_AUC=hit1@metrics$AUC,
          validation_AUC=hit2@metrics$AUC,
          testing_AUC=hit3@metrics$AUC))
  message(paste("Hyperparameters of the best model"))
  print(c(Number_of_trees=rf_model@allparameters$ntrees,
          Mtries=rf_model@allparameters$mtries,
          MaxDepth=rf_model@allparameters$max_depth))
  #h2o.saveModel(rf_model)
  return(rf_model)
}

