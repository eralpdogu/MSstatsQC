#' A data processing function
#'
#' @param data Comma-separated (*.csv), QC file format.
#'  It should contain a Precursor column and the metrics columns.
#' @export
#' @import dplyr
#' @import RecordLinkage
#' @examples
#' # The data is "S9Site54" which is defined in the package.
#' data <- DataProcess(S9Site54)

#########################################################################################################

DataProcess <- function(data = NULL){
  if(is.null(data))
    return()
  ## save process output in each step #### creating a log file ########### from Meena's code
  allfiles <- list.files()

  num <- 0
  filenaming <- "./log/msstatsqc"
  finalfile <- "msstatsqc.log"

  while(is.element(finalfile,allfiles)) {
    num <- num+1
    finalfile <- paste(paste(filenaming,num,sep="-"),".log",sep="")
  }

  #session <- sessionInfo()
  #sink("./log/sessionInfo.txt")
  #print(session)
  #sink()

  #processout <- as.matrix(read.table("./log/sessionInfo.txt", header=TRUE, sep="\t"))
  #write.table(processout, file=finalfile, row.names=FALSE)

  #processout <- rbind(processout, as.matrix(c(" "," ","MSstatsqc - dataProcess function"," "),ncol=1))

  data <- input.sanity.check(data, finalfile)

  return(data)
}
