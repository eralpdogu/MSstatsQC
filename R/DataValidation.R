
#########################################################################################
# here we put a selection of most column names that users use.
#The first element of each vector should be the best name that
# we suggest users to use and  which our code is based on.
#for example "Retention Time" and "Full Width at Half Maximum" which are the first element
# of each vector in the list, are our suggestion so we wrote them in the fisrt place.
best_colnames <- list(
  c("AcquiredTime","Acquired.Time","time"),
  #c("Retention Time","BestRetentionTime" ,"Best.RT","best retention time",
  #  "retention time","rt","best ret time","intensity"),
  #c("Full Width at Half Maximum","MaxFWHM","fwhm","max.fwhm", "Max FWHM"),
  #c("Total Peak Area","Total Area","TotalArea","total area","TA","T.Area"),
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
  # This function removes any existing punctuation in your sentence or word
  #and transfer it to space.
  # Ex1: punc_remove(Best.RT) --> Best RT     #Ex2: punc_remove(Best_RT) --> Best RT
  return(gsub("[[:punct:]///' ]", " ", x))
}
#### clearString function ###############################################################################################
clearString <- function(x){
  # This function, gets a word or setence, Splits it (if it is a camelCase),
  #removes any existing punctuations, and transfer
  # all Upper Case letters to lower case letters.
  # Ex: clearString("myName_isSara.Taheri") --> my name is sara taheri
  return(tolower(punc_remove(camelCaseSplit(x))))
}
#### guessColumnName function ###########################################################################################

# This function receives the data and check the column names of data and changes
#the column names if it is not the
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
#################################################################################################
input.sanity.check <- function(prodata, finalfile) {

  error_message <- ""

  # get the column names and change them to the column names that we want
  #(For ecample we want Retention Time but a user might use RT, this function
  #auotomatically change RT to Retention Time)
  colnames(prodata) <- unlist(lapply(colnames(prodata), function(x)guessColumnName(x)))


  ############## conditions ##############
  # check that the data includes all the requiered columns and if not tell user what column is missing
  required_column_names <- c("Precursor","Annotations")
  if(!("Annotations" %in% colnames(prodata))) {
    prodata[,"Annotations"] <- NA
  }
  provided_column_names <- colnames(prodata)
  if(!all(required_column_names %in% provided_column_names)) {
    missedInput <- which(!(required_column_names %in% provided_column_names))
    error_message <- paste("ERROR : The required input(inputs) : ",
                           paste(required_column_names[missedInput], collapse = ", "),
                           " is(are) not provided in data set. Please add it to your
                           data and try again.\n\n")
  }

  # check that all columns other than Precursor and Acquired Time and
  #Annotations are numeric.
  AfterannoColNum <- (which(colnames(prodata)=="Annotations")) + 1

  for(i in  AfterannoColNum:ncol(prodata)) {
    if(is.numeric(prodata[,i]) == FALSE) {
      error_message <- paste(error_message, "All the values of",
                             colnames(prodata)[i], "should be numeric and positive.\n\n")
    }
  }

  if(error_message != "") {
    return(paste(error_message))
  }
  # for custom metrics we are checking them to be numeric in QCMetrics in
  # "find_custom_metrics" function and only accepting numeric columns after Annotation

  # if there is any missing value in data replace it with NA
  prodata[prodata==""] <- NA

  # Define peak assymetry
  if("MinStartTime" %in% provided_column_names && "MaxEndTime" %in% provided_column_names) {
    peakAss <- 2*prodata$MinStartTime/(prodata$MaxEndTime+prodata$MinStartTime)
    # locate a new column named "Peak Assymetry" right after the column named "Annotation"
    prodata[,"Peak Assymetry"] <- peakAss
  }
  print("Your data is ready to go!")
  return(prodata)
}
