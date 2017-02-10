TestSubsetIdx = function(){
  for(i in 1:10){
    v = SubsetIdx(100,i)
    stopifnot(abs(length(v)-i) <= 2)
    stopifnot(!is.unsorted(v))
  }
}

TestSubsetRoc = function(){
  n = 1000
  x = rnorm(n)
  obs = x + rnorm(n, sd=0.05)
  labels = sign(x)
  roc = ComputeAUC(obs, labels, computeROC=TRUE)$roc
  roc.sub = SubsetRoc(roc, 10)
  stopifnot(roc.sub$tpr %in% roc$tpr)
  stopifnot(roc.sub$fpr %in% roc$fpr)
  stopifnot(!is.unsorted(roc.sub$tpr))
  stopifnot(!is.unsorted(roc.sub$fpr))
}
