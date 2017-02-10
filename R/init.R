
library(affy)
library(devtools)
library(gridExtra)
library(ggplot2)
library(ggrepel)
library(gplots)
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

merge = base::merge
options(error=recover)

### enables list output of function
list = structure(NA,class="result")
"[=.result" = function(x,...,value) {
  args = as.list(match.call())
  args = args[-c(1:2,length(args))]
  length(value) = length(args)
  for(i in seq(along=args)) {
    a = args[[i]]
    if(!missing(a)) eval.parent(substitute(a = v,list(a=a,v=value[[i]])))
  }
  x
}

source('R/src/Utils.R')
source('R/test/test_Utils.R')

source('R/src/DataProc.R')
source('R/test/test_DataProc.R')

source('R/src/TensorROC.R')
source('R/test/test_TensorROC.R')

source('R/src/DefineTensor.R')
source('R/test/test_DefineTensor.R')

source('R/src/GSEA.R')
source('R/test/test_GSEA.R')

source('R/src/Plot.R')
source('R/src/CallMatlab.R')
source('R/src/EvaluateTensor.R')

if(!exists('testAll') || testAll){
  source('R/test/test_all.R')
}

