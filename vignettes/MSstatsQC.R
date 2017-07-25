## ----eval=TRUE-----------------------------------------------------------
source("https://bioconductor.org/biocLite.R")
biocLite("MSstatsQC")

## ---- eval=TRUE----------------------------------------------------------
#A typical multi peptide and multi metric system suitability dataset
#This dataset was generated during CPTAC Study 9.1 at Site 54
library(MSstatsQC)
data <- MSstatsQC::S9Site54

