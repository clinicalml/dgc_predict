
# choose indices between 1 and n_long so that I have approximately n_short
# indices evenly spaced and including the first and last points
SubsetIdx = function(n_long, n_short){
  n_skip = floor(n_long / n_short)
  idx = seq(1, n_long, n_skip)
  if( idx[length(idx)] != n_long){
    idx = c(idx, n_long)
  }
  return(idx)
}

SubsetRoc = function(roc.in, n_short){
  if(length(roc.in) > 2*n_short){
    n_long = length(roc.in$tpr)
    idx = SubsetIdx(n_long, n_short)
    roc.out = list()
    roc.out$tpr = roc.in$tpr[idx]
    roc.out$fpr = roc.in$fpr[idx]
  }else{
    roc.out = roc.in
  }
  return(roc.out)
}


