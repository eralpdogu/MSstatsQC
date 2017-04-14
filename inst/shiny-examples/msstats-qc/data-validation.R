COL.BEST.RET <- "Retention Time"
COL.FWHM <- "Full Width at Half Maximum"
COL.TOTAL.AREA <- "Total Peak Area"
COL.PEAK.ASS <- "Peak Assymetry"
#########################################################################################
# here we put a selection of most column names that users use. The first element of each vector should be the best name that
# we suggest users to use and  which our code is based on. for example "Retention Time" and "Full Width at Half Maximum" which are the first element
# of each vector in the list, are our suggestion so we wrote them in the fisrt place.
best_colnames <- list(
  c("AcquiredTime","Acquired.Time","time","creation date"),
  c("Retention Time","BestRetentionTime" ,"Best.RT","best retention time", "retention time","rt","best ret time","intensity","Best RT"),
  c("Full Width at Half Maximum","MaxFWHM","fwhm","max.fwhm", "Max FWHM"),
  c("Total Peak Area","Total Area","TotalArea","total area","TA","T.Area"),
  c("MinStartTime","min start time","Min Start Time"),
  c("MaxEndTime", "max end time","Max End Time"),
  c("Precursor","PeptideSequence"),
  c("Annotations","anotations")
)
#### camelCaseSplit function ##############################################################################################
camelCaseSplit <- function(x) {
  # This function get a camelCase word and splits it.
  # Ex : camelCaseSplit("myComputerIsHere") ---> my Computer Is Here
  return(gsub("([a-z])([A-Z])", "\\1 \\L\\2", x, perl = TRUE))
}
#### punc_remove function #################################################################################################
punc_remove <- function(x){
  # This function removes any existing punctuation in your sentence or word and transfer it to space.
  # Ex1: punc_remove(Best.RT) --> Best RT     #Ex2: punc_remove(Best_RT) --> Best RT
  return(gsub("[[:punct:]///' ]", " ", x))
}
#### clearString function ###############################################################################################
clearString <- function(x){
  # This function, gets a word or setence, Splits it (if it is a camelCase), removes any existing punctuations, and transfer
  # all Upper Case letters to lower case letters.
  # Ex: clearString("myName_isSara.Taheri") --> my name is sara taheri
  return(tolower(punc_remove(camelCaseSplit(x))))
}
#### guessColumnName function ###########################################################################################

# This function receives the data and check the column names of data and changes the column names if it is not the
# same names as our suggested sample data to fit our suggested sample data
guessColumnName <- function(x){

a <- clearString(x)

  max_index <- 0
  max <- -1
  for(i in 1:length(best_colnames)){
    col <- best_colnames[[i]]
    for(j in 1:length(col)){
      sim <- levenshteinSim(a,clearString(col[j]))
      if(sim > max){
        max <- sim
        max_index <- i
      }
    }
  }
  if (max > 0.6) {
    return(best_colnames[[max_index]][1])
  }
  else {
    return(x)
  }
}

input.sanity.check <- function(prodata, processout, finalfile) {
  error_message <- ""

  # get the column names and change them to the column names that we want (For ecample we want Retention Time but a user might use RT, this function auotomatically change RT to Retention Time)
  colnames(prodata) <- unlist(lapply(colnames(prodata), function(x)guessColumnName(x)))


  ### conditions
  # check that the data includes all the requiered columns and if not tell user what column is missing
  # required_column_names <- c("Precursor","Retention Time","Full Width at Half Maximum","Total Peak Area","MinStartTime"
  #                            ,"MaxEndTime")
  required_column_names <- c("Precursor","Annotations")
  if(!("Annotations" %in% colnames(prodata))) {
    prodata[,"Annotations"] <- NA
  }
  provided_column_names <- colnames(prodata)
  # if(!all(required_column_names %in% provided_column_names)) {
  #   missedInput <- which(!(required_column_names %in% provided_column_names))
  #   error_message <- paste("ERROR : The required input(inputs) : ",
  #                          paste(required_column_names[missedInput], collapse = ", "),
  #                          " is(are) not provided in data set. Please add it to your data and try again.\n\n")
  # }

  # check that all columns other than Precursor and Acquired Time and Annotations are numeric.
  AfterannoColNum <- (which(colnames(prodata)=="Annotations")) + 1

  for(i in  AfterannoColNum:ncol(prodata)) {
    if(is.numeric(prodata[,i]) == FALSE) {
      error_message <- paste(error_message, "All the values of", colnames(prodata)[i], "should be numeric and positive.\n\n")
    }
  }
   # if(!is.numeric(prodata[,COL.BEST.RET]) || !all(prodata[,COL.BEST.RET] > 0)) {
   #   error_message <- paste(error_message, "All the values of Retention Time should be numeric and positive.\n\n")
   #   #return(error_message)
   # }


   # if(!is.numeric(prodata[,COL.FWHM]) || !all(prodata[,COL.FWHM] > 0)) {
   #   error_message <- paste(error_message,"All the values of Full Width at Half Maximum should be numeric and positive.\n\n")
   #   #return(error_message)
   # }

  # if(!is.numeric(prodata[,COL.TOTAL.AREA]) || !all(prodata[,COL.TOTAL.AREA] > 0)) {
  #   error_message <- paste(error_message,"All the values of Total Peak Area should be numeric and positive.\n\n")
  #   #return(error_message)
  # }


  # if(!is.numeric(prodata$MaxEndTime) || !all(prodata$MaxEndTime > 0)) {
  #   error_message <- paste(error_message,"All the values of Max End Time should be numeric and positive.\n\n")
  #   #return(error_message)
  # }

  # if(!is.numeric(prodata$MinStartTime) || !all(prodata$MinStartTime > 0)) {
  #   error_message <- paste(error_message,"All the values of Min End Time should be numeric and positive. \n\n")
  #   #return(error_message)
  # }

  if(error_message != "") {
    return(paste(error_message, "Please check the values to make sure all the inputs are numeric and positive and then try again."))
  }
  # for custom metrics we are checking them to be numeric in QCMetrics in "find_custom_metrics" function and only accepting numeric columns after Annotation

  # if there is any missing value in data replace it with NA
  prodata[prodata==""] <- NA
  # some times numeric values of some users are like 333,222 which is not acceptable and we convert it to 333222 by replacing "," to ""
  # prodata[,"Full Width at Half Maximum"] <- as.numeric(gsub(",","",prodata[,"Full Width at Half Maximum"]))
  # prodata[,"Total Peak Area"] <- as.numeric(gsub(",","",prodata[,"Total Peak Area"]))
  # prodata[,"Retention Time"] <- as.numeric(gsub(",","",prodata[,"Retention Time"]))
  # prodata$MaxEndTime <- as.numeric(gsub(",","",prodata$MaxEndTime))
  # prodata$MinStartTime <- as.numeric(gsub(",","",prodata$MinStartTime))


  # some data migh have annotation column, some might not have. If it doesn't, we create an empty "Annotation" column at the very end column of the data


  # Define peak assymetry
  if("MinStartTime" %in% provided_column_names && "MaxEndTime" %in% provided_column_names) {
    peakAss <- 2*prodata$MinStartTime/(prodata$MaxEndTime+prodata$MinStartTime)
    # locate a new column named "Peak Assymetry" right after the column named "Annotation"
    #prodata.first <- prodata[,1:which(colnames(prodata)=="Annotations")]
    #prodata.first[,"Peak Assymetry"]<- peakAss
    #prodata <- cbind(prodata.first, prodata[,(which(colnames(prodata)=="MaxEndTime")+1):ncol(prodata), drop = FALSE])
    prodata[,"Peak Assymetry"] <- peakAss
  }

  return(prodata)

}

### Input_checking function #########################################################################################
input_checking <- function(data){

  data <- input.sanity.check(data, processout, finalfile)


  return(data)
}
