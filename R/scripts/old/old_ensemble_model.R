

# LocalEnsEstimate_DrugCell = function(T_knn, T_tensor, lambda=0.5){
#   # Get max number of cell types available per drug
#   nPerDrug = NumSigs(T_knn, 'drug')
#   Nc = max(nPerDrug)
#   
#   # Get max number of drugs available per cell type
#   nPerCell = NumSigs(T_knn, 'cell')
#   Nd = max(nPerCell)
#   
#   ens = array(data=NaN, dim=dim(T_knn), dimnames=dimnames(T_knn))
#   lambdas = c()
#   cnt = 1
#   
#   for(i in 1:dim(T_knn)[1]){
#     print(i)
#     f1 = nPerDrug[i] / Nc
#     for(k in 1:dim(T_knn)[3]){
#       if(!is.na(T_knn[i,1,k])){
#         f2 = nPerCell[k] / Nd
#         lambda_ik = (lambda * f1 + f2) / (lambda + 1)
#         ens[i,,k] = lambda_ik * T_tensor[i,,k] + (1-lambda_ik) * T_knn[i,,k]
#         lambdas[cnt] = lambda_ik
#         cnt = cnt + 1
#       }
#     }
#   }
#   
#   return(list(ens=ens, lambdas=lambdas))
# }

# idx = c(1:50, 2081:2130)
# sm = list(dnpp=tensors$cv$dnpp[idx,,1:10], tensor=tensors$cv$tensor[idx,,1:10], meas=tensors$meas[idx,,1:10])
# sm_out = LocalEnsEst_Drug(sm_tensors$dnpp, sm_tensors$tensor, lambda=0.5)
# #test_local = sapply(seq(0, 1, 0.05), function(lambda){print(lambda); tmp = LocalEnsEstimate(sm$dnpp, sm$tensor, lambda); ComputePCT(sm$meas, tmp$ens)})
# i=1
# tmp = list()
# test_local = list()
# for(lambda in alphas){
#   print(lambda)
#   tmp[[i]] = LocalEnsEst_Drug(sm$dnpp, sm$tensor, lambda)
#   test_local[[i]] = ComputePCT(sm$meas, tmp[[i]]$ens)
#   i = i+1
# }
