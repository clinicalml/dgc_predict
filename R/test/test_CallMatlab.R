TestStartMatlab = function(){
  # See TestMatlab below, this essentially tests whether Matlab has been started appropriately
}

TestMatlab = function(mat=get0('matlab')){
  working = FALSE
  if(!is.null(mat) && isOpen(mat)){
    evaluate(mat, 'a=2+2;')
    a = getVariable(mat, 'a')$a
    if(a == 4){
      working = TRUE
    }
  }
  return(working)
}

TestCompleteTensor= function(mat = get0('matlab')){
  matlab = StartMatlab()
  tensorFull = MakeTensor(NA_frac=0, fill=5)
  tensorNA = tensorFull
  tensorNA[5,,1] = NA
  tensorComplete = CompleteTensor(matlab, tensorNA, method='mean')
  stopifnot(identical(tensorComplete, tensorFull))
  
  tensorComplete = CompleteTensor(matlab, tensorNA, method='mean2')
  stopifnot(identical(tensorComplete, tensorFull))
  
  tensorComplete = CompleteTensor(matlab, tensorNA, method='knnd')
  stopifnot(identical(tensorComplete, tensorFull))
  
  tensorComplete = CompleteTensor(matlab, tensorNA, method='fa_lrtc')
  stopifnot(Norm2(tensorFull[5,,1]- tensorComplete[5,,1]) < 0.1)
  
  close(matlab)
  rm(matlab)
}


TestCrossValidateTensor = function(mat = get0('matlab')){
  matlab = StartMatlab()
  tensor = NormSigs(MakeTensor(nDrugs=50, nCells=10, nGenes=50, NA_frac=0, fill=5, removeSig_frac=0.2, noise=0.1))
  out = CrossValidateTensor(matlab, tensor, methods=c('mean'))
  close(matlab)
  rm(matlab)
  stopifnot(abs(ComputePCT(out$tensors, tensor) - out$PCT) < 1e-5)
}
