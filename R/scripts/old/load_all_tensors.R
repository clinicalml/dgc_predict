
library(rhdf5)

# Read all tensor data from large tensor experiments

tensors = list()

print('Loading data tensor..')
tensors$meas = LoadTensorMat(DataDir(sprintf('tensors/%s.mat', tsize)))$tensor
names(dimnames(tensors$meas)) = c('drug','gene','cell')

print('Loading cross-validated tensors..')
file = ResultsDir(sprintf('%s/%s_tensor_results.mat', tsize, tsize))
tensors$cv = list(mean=h5read(file, '#refs#/b'), mean2=h5read(file,'#refs#/c'),
                  dnpp=h5read(file,'#refs#/d'), tensor=h5read(file,'#refs#/e'))
tensors$cv = lapply(tensors$cv, function(tensor){dimnames(tensor) = dimnames(tensors$meas); return(tensor)})

if(loadMergeAndPred){
  print('Loading completed tensors and then generating merged and predicted tensors...')
  for(method in c('mean', 'mean2','dnpp', 'tensor')){
    print(method)
    
    print('..completed tensor')
    if(tsize == 'large'){
      tcomp = h5read(ResultsDir(sprintf('large/hdf5/%s_final_hdf5.mat', method)),'T')
    }else{
      tcomp = LoadTensorMat(DataDir(sprintf('results/tsize/%s/%s_final.mat', tsize, method)))$tensor
    }
    
    dimnames(tcomp) = dimnames(tensors$meas)
    tcomp = NormSigs(tcomp)
    
    print('..merged tensor')
    mergeT = tensors$cv[[method]]
    mergeT[is.na(mergeT)] = tcomp[is.na(mergeT)]
    tensors$merge[[method]] = mergeT
    
    print('..predicted tensor')
    predT = tcomp
    predT[!is.na(tensors$cv[[method]])] = NA
    tensors$pred[[method]] = predT
    stopifnot(NumSigs(tensors$pred[[method]]) + NumSigs(tensors$cv[[method]]) == NumSigs(tcomp))
    
    rm(predT, tcomp, mergeT)
  }
}
