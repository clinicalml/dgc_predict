

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