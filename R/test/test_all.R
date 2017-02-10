source('R/test/test_helpers.R')
debug=0

function.names = c(GetFunctionNames('R/src/Utils.R'), 
                    GetFunctionNames('R/src/DataProc.R'),
                    GetFunctionNames('R/src/TensorROC.R'),
                    GetFunctionNames('R/src/DefineTensor.R'))

for (i in 1:length(function.names)){
  func = paste('Test', function.names[i], '()', sep='')
  cat(paste(func, '...', sep=''))
  eval(parse(text=func))
  cat(' passed\n')
}
