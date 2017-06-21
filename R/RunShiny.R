#' A function to launch MSstatsQC Shiny app.

#' @param data Comma-separated (*.csv), QC file format. It should contain a Precursor column
#'  and the metrics columns.
#' @keywords shiny
#' @export
#' @import shiny
#' @examples
#' sampleData <- DataProcess(S9Site54)
#' RunShiny(data = sampleData)

RunShiny <- function(data = NULL) {
  appDir <- system.file(package = "MSstatsQC")
  source(paste0(appDir,"/shiny-examples/msstats-qc/app.R"))
  runner(data) # to be implemented
}
