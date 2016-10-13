
library(plyr)
library(R.matlab)

options(error=recover)
 
### enables list output of function
list <- structure(NA,class="result")
"[<-.result" <- function(x,...,value) {
  args <- as.list(match.call())
  args <- args[-c(1:2,length(args))]
  length(value) <- length(args)
  for(i in seq(along=args)) {
    a <- args[[i]]
    if(!missing(a)) eval.parent(substitute(a <- v,list(a=a,v=value[[i]])))
  }
  x
}

source('R/utils/Utils.R')
source('R/utils/test_Utils.R')

source('R/dataprocessing/DataProc.R')
source('R/dataprocessing/test_DataProc.R')

source('R/analyze_tensor_results/TensorROC.R')
source('R/analyze_tensor_results/test_TensorROC.R')

source('R/drug_repurp/DrugRepurp.R')
source('R/drug_repurp/test_DrugRepurp.R')

source('R/dataprocessing/DefineTensor.R')
source('R/dataprocessing/test_DefineTensor.R')

source('R/gsea/GSEA.R')
source('R/gsea/test_GSEA.R')

source('R/plot/Plot.R')
source('R/call_matlab/CallMatlab.R')
source('R/analyze_tensor_results/EvaluateTensor.R')

if(!exists('testAll') || testAll){
  source('R/test/test_all.R')
}

