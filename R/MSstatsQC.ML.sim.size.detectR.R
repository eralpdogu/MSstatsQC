#' A function to train random forest classifiers for QC data
#'
#' @param guide.set comma-separated (.csv), metric file. It should contain a "Precursor" column and the metrics columns. It should also include "Annotations" for each run.
#' @param sim.start enter min simulation size.
#' @param sim.end enter max simulation size.
#' @export
#' @import h2o
#' @import ggplot2
#' @import gridExtra
#' @import FrF2
#' @import car
#' @import reshape2
#' @import dplyr
#' @import plyr
#' @return a plot for sim.size vs performance
#' @keywords machine learning, simulation size
#' @examples
#' # First process the data to make sure it's ready to use
#' S9Site54.dataML<-DataProcess(MSstatsQC::S9Site54[,])
#' colnames(S9Site54.dataML)[1]<-c("idfile")
#' colnames(S9Site54.dataML)[2]<-c("peptide")
#' S9Site54.dataML$peptide <-as.factor(S9Site54.dataML$peptide)
#' S9Site54.dataML$idfile<-as.numeric(S9Site54.dataML$idfile)
#' S9Site54.dataML <- within(S9Site54.dataML, rm(Annotations, missing))
#' guide.set <- filter(S9Site54.dataML, idfile<=20)
#' MSstatsQC.ML.sim.size.detectR(guide.set, sim.start=10, sim.end=2500)


MSstatsQC.ML.sim.size.detectR<-function(guide.set, sim.start, sim.end){

  sequence<-seq(sim.start, sim.end, 20)
  results<-matrix(NA, 100000, 4)
  for(i in sequence){
    rf_model<-MSstatsQC.ML.trainR(guide.set, i, guide.set.annotations = NULL)
    cf<- data.frame(h2o.confusionMatrix(rf_model),stringsAsFactors = F)
    sens<-cf[1,1]/(cf[1,1]+cf[2,1])
    err1<-cf[1,3]
    err3<-cf[3,3]
    results[i,]<-cbind(i,sens,err1,err3)
  }
  results<-as.data.frame(results[complete.cases(results),])
  colnames(results)<-c("Simulation.size", "Accuracy", "False positive rate", "False negative rate")

  results_melt <- melt(results,id.vars ="Simulation.size",
                       measure.vars=c("Accuracy", "False positive rate", "False negative rate"))
  ggplot(results_melt, aes(Simulation.size, value)) +
    geom_point()+
    geom_smooth(aes(color="red"), show.legend = FALSE)+
    ylab("Probability")+
    xlab("Simulation size")+
    facet_wrap(~variable,scales = "free")+
    theme(text = element_text(size=12),
          axis.text.x = element_text(angle=45, hjust=1))
}

