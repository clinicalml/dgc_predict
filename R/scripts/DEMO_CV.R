tensors = list()
tensors$meas = LoadTensorMat(DataDir('tensors/T50_1.RData'))

matlab = StartMatlab()

for(method in c('mean', 'mean2', 'knn', 'fa_lrtc')){
  if(method == 'fa_lrtc'){
    mthd = 'tensor'
  }else{mthd = method}
  tensors$cv[[mthd]] = CrossValidateTensor(matlab, tensors$meas, method, dataset, 
                                           nFolds = NumSigs(tensors$meas), 
                                           maxFolds = NumSigs(tensors$meas))
  dimnames(tensors$cv[[mthd]]) = dimnames(tensors$meas)
}
close(matlab)
rm(matlab)