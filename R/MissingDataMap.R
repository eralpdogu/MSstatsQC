#' A function to summarize missing values
#'
#' @param data Processed data
#' @return A plot of missing values.
#' @export
#' @import ggplot2
#' @examples
#' # The data is "S9Site54" which is defined in the package.
#' #data <- DataProcess(S9Site54)
#' #MissingDataMap(data)

MissingDataMap <- function(data){

  if (sum(data$missing)==0) {
    print("No missing values!")
    }
  else {
  ggplot(data, aes(x=AcquiredTime, y=Precursor))+geom_raster(aes(fill = missing))+
    xlab("Time")+
    ylab("Peptide")+
    scale_fill_continuous("# of missing metrics", low = "lightblue", high = "white")+
    theme_bw()+
    theme(legend.position="bottom")+
    theme( panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           text = element_text(size=10),
           axis.text.x = element_text(angle=45, hjust=1))
  }
}



