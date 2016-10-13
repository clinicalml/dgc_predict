source('R/test/test_helpers.R')
debug=0

function.names <- c(GetFunctionNames('R/utils/Utils.R'), 
                    GetFunctionNames('R/dataprocessing/DataProc.R'),
                    GetFunctionNames('R/analyze_tensor_results/TensorROC.R'),
                    GetFunctionNames('R/dataprocessing/DefineTensor.R'),
                    GetFunctionNames('R/drug_repurp/DrugRepurp.R'))

for (i in 1:length(function.names)){
  func <- paste('Test', function.names[i], '()', sep='')
  cat(paste(func, '...', sep=''))
  eval(parse(text=func))
  cat(' passed\n')
}
