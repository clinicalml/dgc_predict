StartMatlab = function(mat = get0('matlab')){
  if(is.null(mat)){
    library(R.matlab)
    Matlab$startServer(matlab=GetConfig('MTLABPATH')) 
    mat = Matlab()
  }
  if(!isOpen(mat)){
    isOpen = open(mat)
    if (!isOpen) throw('MATLAB server is not running: waited 30 seconds.')
  }
  return(mat)
}

CompleteTensor = function(matlab, tensor, method){
  setVariable(matlab, method=method, T=tensor)
  evaluate(matlab, 'args = GetArgs(method, [], [], size(T));')
  evaluate(matlab, 'out = CompleteTensor(T, method, args);')
  out = getVariable(matlab, 'out')$out
  dimnames(out) = dimnames(tensor)
  return(out)
}

CrossValidateTensor = function(matlab, tensor, methods=c('mean','mean2', 'knnd', 'fa_lrtc'), exp_name='test', nFolds=10, maxFolds=10, saveFile=FALSE){
  setVariable(matlab, methods=methods, tensor=tensor, exp_name=exp_name, nFolds=nFolds, maxFolds=maxFolds)
  evaluate(matlab, '[cvTensors, PCT, PCTf] = TensorCV4(methods, tensor, exp_name, nFolds, maxFolds);')
  out = getVariable(matlab, 'cvTensors')$cvTensors
  PCT = as.vector(getVariable(matlab, 'PCT')$PCT)
  PCTf = getVariable(matlab, 'PCTf')$PCTf
  n = length(methods)
  stopifnot(length(out) == n)
  if(n > 1){
    outR = lapply(1:n, function(i){A = out[[i]][[1]]; dimnames(A) = dimnames(tensor); return(A)})
    names(outR) = methods
    colnames(PCTf) = methods
  }else{
    outR = out[[1]][[1]]
    dimnames(outR) = dimnames(tensor)
    PCTf = as.vector(PCTf)
  }
  names(PCT) = methods
  return(list(tensors=outR, PCT=PCT, PCTf=PCTf))
}