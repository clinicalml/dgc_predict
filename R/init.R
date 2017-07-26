
library(affy)
library(data.table)
library(devtools)
library(fields)
library(gridExtra)
library(ggplot2)
library(ggrepel)
library(gplots)
library(GSEABase)
library(GO.db)
library(HTSanalyzeR)
library(KEGG.db)
library(matrixStats)
library(org.Hs.eg.db)
library(plotrix)
library(plyr)
library(RColorBrewer)
library(reshape)
library(reshape2)
library(ROCR)
library(R.matlab)
library(R.utils)
library(WriteXLS)

# These are from GSEA.R. Not sure which ones I need.
# library(cellHTS2)
# library(VennDiagram)
# library(RCurl)
# library(snow)

merge = base::merge
options(error=recover)

### enables list output of function
list <- structure(NA,class="result")
"[<-.result" <- function(x,...,value) {
  args = as.list(match.call())
  args = args[-c(1:2,length(args))]
  length(value) = length(args)
  for(i in seq(along=args)) {
    a = args[[i]]
    if(!missing(a)) eval.parent(substitute(a <- v,list(a=a,v=value[[i]])))
  }
  x
}

source('R/src/Utils.R')
source('R/test/test_Utils.R')

source('R/src/DataProc.R')
source('R/test/test_DataProc.R')

source('R/src/DefineTensor.R')
source('R/test/test_DefineTensor.R')

source('R/src/EvaluateTensor.R')
source('R/test/test_EvaluateTensor.R')

source('R/src/CallMatlab.R')
source('R/test/test_CallMatlab.R')

source('R/src/GSEA.R')
source('R/test/test_GSEA.R')

source('R/src/Ensemble.R')
#source('R/test/test_Ensemble.R')

source('R/src/Plot.R')

if(!exists('testAll') || testAll){
  source('R/test/test_all.R')
}


