tensors = LoadTensors(tsize='small', print=FALSE)
T_knn = tensors$cv$dnpp
T_tensor = tensors$cv$tensor
manycell = LoadTensorMat(DataDir('tensors/manycell.mat'))$tensor

TestEnsembleGlobal = function(){
  CheckTensors()
  
  out = EnsembleGlobal(T_knn=T_knn, T_tensor=T_tensor, lambda=0)
  stopifnot(identical(out, T_tensor))
  
  out = EnsembleGlobal(T_knn=T_knn, T_tensor=T_tensor, lambda=1)
  stopifnot(identical(out, T_knn))
  
  out = EnsembleGlobal(T_knn=T_knn, T_tensor=T_tensor, lambda=0.5)
  stopifnot(out[4, 20, 5] == mean(c(T_knn[4, 20, 5], T_tensor[4, 20, 5])))
}

TestEnsembleLocal = function(){
  CheckTensors()
  
  out0 = EnsembleLocal(T_knn=T_knn, T_tensor=T_tensor, lambda=0)
  stopifnot(identical(out0$ens, T_knn))
  
  lambdas = c(0.5, 1, 2, 5)
  out = lapply(lambdas, function(lam) EnsembleLocal(T_knn=T_knn, T_tensor=T_tensor, lambda=lam))
  pct = sapply(out, function(x) ComputePCT(T_knn, x$ens))

  stopifnot(all(sort(pct, decreasing=TRUE) == pct))
  
  idx = which.max(NumSigs(T_knn, 'cell'))
  tmp = lapply(out, function(x) stopifnot(identical(T_knn[,,idx], x$ens[,,idx])))
}

TestEnsembleHeuristic = function(){
  CheckTensors()
  
  # For small tensor, the entire thing should be T_tensor
  out = EnsembleHeuristic(T_knn=T_knn, T_tensor=T_tensor, manycell=manycell)
  stopifnot(identical(T_tensor, out))
  
  # # But this is really meant to be applied to the large tensor (commented out because it takes a long time)
  # tlg = LoadTensors(tsize='large', print=TRUE)
  # out = EnsembleHeuristic(T_knn=tlg$cv$dnpp, T_tensor=tlg$cv$tensor, manycell=manycell)
  # mc = NumSigs(manycell)
  # n_knn = length(which(tlg$cv$dnpp[,1,] == out[,1,])) 
  # n_tensor = length(which(tlg$cv$tensor[,1,] == out[,1,]))
  # stopifnot(abs(n_tensor - mc) < 10)
  # stopifnot(abs(n_knn - (NumSigs(tlg$cv$dnpp) - mc)) < 10)
}

CheckTensors = function(){
  stopifnot(dim(T_knn) == c(300, 978, 15) && NumSigs(T_knn) == 3210)
  stopifnot(dim(T_tensor) == c(300, 978, 15) && NumSigs(T_tensor == 3210))
}



