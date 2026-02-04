#' A data processing function for removing missing values
#'
#' @param data Comma-separated (*.csv), QC file format. It should contain a Precursor column and the metrics columns.
#' @return A data frame that processes using \code{input.sanity.check} function.
#' @export
#' @importFrom stats complete.cases
#' @examples
#' # The data is "S9Site54" which is defined in the package.
#' data <- RemoveMissing(S9Site54)
RemoveMissing <- function(data = NULL) {
    data <- data[complete.cases(data), ] # work with complete cases
    return(data)
}
