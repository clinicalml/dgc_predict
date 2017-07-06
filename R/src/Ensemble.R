
EnsembleGlobal = function(T_knn, T_tensor, lambda=0.5){
  stopifnot(lambda <= 1 && lambda >= 0)
  return(lambda*T_knn + (1-lambda)*T_tensor)
}

EnsembleLocal = function(T_knn, T_tensor, lambda=0.5){
  # Get fraction of drugs available per cell type
  nPerCell = NumSigs(T_knn, 'cell')
  Nd = max(nPerCell)
  drugFracPerCell = nPerCell / Nd
  
  ens = array(data=NaN, dim=dim(T_knn), dimnames=dimnames(T_knn))
  lambdas = c()
  cnt = 1
  
  for(c in 1:dim(T_knn)[3]){ # for each cell type
    lambda_c = drugFracPerCell[c] ^ lambda
    for(d in 1:dim(T_knn)[1]){ # for each drug
      if(!is.na(T_knn[d,1,c])){
        ens[d,,c] = lambda_c * T_knn[d,,c] + (1-lambda_c) * T_tensor[d,,c]
        lambdas[cnt] = lambda_c
        cnt = cnt + 1
      }
    }
  }
  
  stopifnot(NumSigs(ens) == NumSigs(T_tensor))
  
  return(list(ens=ens, lambdas=lambdas))
}

EnsembleHeuristic = function(T_knn, T_tensor, manycell){
  drugs = dimnames(manycell)[[1]]
  cells = dimnames(manycell)[[3]]
  ens = T_knn
  for(drug in dimnames(T_knn)[[1]]){
    for(cell in dimnames(T_knn)[[3]]){
      if(drug %in% drugs & cell %in% cells){
        ens[drug,,cell] = T_tensor[drug,,cell]  
      }
    }
  }
  return(ens)
}
