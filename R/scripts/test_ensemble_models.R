# Load all tensors (input and cross-validated)
tensors = LoadTensors(tsize='small', print=TRUE)

lambdas = seq(0, 1, 0.05)

# Method 1: Simply estimate a global trade-off parameter, lambda
pct1 = sapply(lambdas, function(lambda){print(lambda); ens = GlobalEnsEstimate(tensors$cv$dnpp, tensors$cv$tensor, lambda); ComputePCT(tensors$meas, ens)})
err1 = sapply(lambdas, function(lambda){print(lambda); ens = GlobalEnsEstimate(tensors$cv$dnpp, tensors$cv$tensor, lambda); ComputeErrorRate(tensors$meas, ens)})

# Method 2: Estimate a local trade-off parameter, lambda
pct2 = sapply(lambdas, function(lambda){print(lambda); tmp = LocalEnsEst_Drug(tensors$cv$dnpp, tensors$cv$tensor, lambda); ComputePCT(tensors$meas, tmp$ens)})





manycell = LoadTensorMat(DataDir('tensors/manycell.mat'))$tensor


idx = c(1:50, 2081:2130)
sm = list(dnpp=tensors$cv$dnpp[idx,,1:10], tensor=tensors$cv$tensor[idx,,1:10], meas=tensors$meas[idx,,1:10])
#sm = list(dnpp=tensors$cv$dnpp[1:100,,1:10], tensor=tensors$cv$tensor[1:100,,1:10], meas=tensors$meas[1:100,,1:10])
sm_out = LocalEnsEst_Drug(sm_tensors$dnpp, sm_tensors$tensor, lambda=0.5)
#test_local = sapply(seq(0, 1, 0.05), function(lambda){print(lambda); tmp = LocalEnsEstimate(sm$dnpp, sm$tensor, lambda); ComputePCT(sm$meas, tmp$ens)})
i=1
tmp = list()
test_local = list()
for(lambda in alphas){
  print(lambda)
  tmp[[i]] = LocalEnsEst_Drug(sm$dnpp, sm$tensor, lambda)
  test_local[[i]] = ComputePCT(sm$meas, tmp[[i]]$ens)
  i = i+1
}

