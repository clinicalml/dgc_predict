
# Run parameters
tsize = 'large'

# Load all tensors (input and cross-validated)
tensors = LoadTensors(tsize=tsize, print=TRUE)


# Method 1: Simply estimate a global trade-off parameter, lambda
param_global = seq(0, 1, 0.05)
global = sapply(param_global, function(lambda){
  print(lambda)
  ens = EnsembleGlobal(T_knn=tensors$cv$dnpp, T_tensor=tensors$cv$tensor, lambda)
  pct = ComputePCT(tensors$meas, ens)
  err = ComputeErrorRate(tensors$meas, ens)
  return(list(pct=pct, err=err))})


# Method 2: Estimate a local trade-off parameter lambda
# The "localness" is that the trade-off parameter is specific to each cell type, depending on how many different drug signatures are available.
param_local = seq(-0.2, 2.0, 0.1)
local = sapply(param_local, function(lambda){
  print(lambda)
  tmp = EnsembleLocal(T_knn=tensors$cv$dnpp, T_tensor=tensors$cv$tensor, lambda)
  pct = ComputePCT(tensors$meas, tmp$ens)
  err = ComputeErrorRate(tensors$meas, tmp$ens)
  return(list(pct=pct, err=err))})


# Method 3: Use a heuristic where we use the tensor approach in the ManyCell region, otherwise use KNN
manycell = LoadTensorMat(DataDir('tensors/manycell.mat'))$tensor
est_heur = EnsembleHeuristic(T_knn=tensors$cv$dnpp, T_tensor=tensors$cv$tensor, manycell=manycell)
heuristic = list(pct=ComputePCT(tensors$meas, est_heur), err=ComputeErrorRate(tensors$meas, est_heur))


# Compare the three approaches:

par(mfrow=c(2,2))
pct_ylim = range(range(global[1,]), range(local[1,]))
err_ylim = range(range(local[2,]), range(global[2,]))
glob_xlim = range(param_global)
loc_xlim = range(param_local)

# Global method, PCT
idx = which.max(global[1,])
plot(param_global, global[1,], xlab='Weight on KNN Estimate', ylab='PCT Accuracy', type='l', lwd=2, xlim=glob_xlim, ylim=pct_ylim,
     main=sprintf('Global ensemble approach: \nAccuracy on %s tensor (optimal = %0.2f)', tsize, param_global[idx])) 
par(new=TRUE)
plot(param_global[idx], global[1,idx], col='red', lwd=2, pch=8, xlim=glob_xlim, ylim=pct_ylim, xlab='', ylab='')
abline(h=heuristic$pct, lwd=2, col='blue')

# Local method, PCT
idx = which.max(local[1,])
plot(param_local, local[1,], xlab='Alpha (exponent on measured drug fraction)', ylab='PCT Accuracy', type='l', lwd=2, xlim=loc_xlim, ylim=pct_ylim,
     main=sprintf('Local ensemble approach: \nAccuracy on %s tensor (optimal = %0.2f)', tsize, param_local[idx]))
par(new=TRUE)
plot(param_local[idx], local[1,idx], col='red', lwd=2, pch=8, xlim=loc_xlim, ylim=pct_ylim, xlab='', ylab='')
abline(h=heuristic$pct, lwd=2, col='blue')

# Global method, Error Rate
idx = which.min(global[2,])
plot(param_global, global[2,], xlab='Weight on KNN Estimate', ylab='Error Rate', type='l', lwd=2, xlim=glob_xlim, ylim=err_ylim,
     main=sprintf('Global ensemble approach: \nErrror rate on %s tensor (optimal = %0.2f)', tsize, param_global[idx])) 
par(new=TRUE)
plot(param_global[idx], global[2,idx], col='red', lwd=2, pch=8, xlim=glob_xlim, ylim=err_ylim, xlab='', ylab='')
abline(h=heuristic$err, lwd=2, col='blue')

# Local method, error rate
idx = which.min(local[2,])
plot(param_local, local[2,], xlab='Alpha (exponent on measured drug fraction)', ylab='Error Rate', type='l', lwd=2, xlim=loc_xlim, ylim=err_ylim,
     main=sprintf('Local ensemble approach: \nError rate on %s tensor (optimal = %0.2f)', tsize, param_local[idx]))
par(new=TRUE)
plot(param_local[idx], local[2,idx], col='red', lwd=2, pch=8, xlim=loc_xlim, ylim=err_ylim, xlab='', ylab='')
abline(h=heuristic$err, lwd=2, col='blue')

