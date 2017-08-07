files = c('Utils', 'DataProc', 'DefineTensor', 'EvaluateTensor', 'CallMatlab', 'GSEA')
functionNames = as.character(unlist(sapply(files, function(file) GetFunctionNames(sprintf('R/src/%s.R', file)))))

for(func in functionNames){
 funcName = paste('Test', func, '()', sep='')
 cat(paste(func, '...', sep=''))
 eval(parse(text=funcName))
 cat(' passed\n')
}
