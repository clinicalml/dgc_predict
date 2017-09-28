TestComputeGCP = function(){
  set.seed(123)
  T1 = MakeTensor()
  T2 = list(test_tensor = MakeTensor())
  GCP = ComputeGCP(T1, T2, plot=FALSE, subset='cv', file=NA, print=FALSE)
  
  # Test: should be looking at correlations between a bunch of random small
  # numbers, and hence GCP should be small
  stopifnot(all(abs(GCP$test_tensor) < 0.2))
}

TestComputePCT = function(){
  set.seed(123)
  T1 = MakeTensor()
  T2 = T1
  
  # Check that even if I remove different random elements from each of T1 and
  # T2, I still get perfect correlation based on remaining entries
  T1 = RemoveRandomEltsFromTensor(T1, NA_frac=0.1)
  T2 = RemoveRandomEltsFromTensor(T2, NA_frac=0.1)

  list[x1, x2] = Tensor2Vec(T1, T2)
  stopifnot(length(x1) < prod(dim(T1))*0.9)
  
  stopifnot(ComputePCT(T1, T2) == 1)
  
  # Check that a constant difference doesn't affect the result
  T2 = T2 - 100
  stopifnot(ComputePCT(T1, T2) == 1)
  
  # But a negative multiplier does
  T2 = -3*T2
  stopifnot(ComputePCT(T1, T2) == -1)
}

TestComputePCT_AllModes = function(){
  set.seed(123)
  T1 = MakeTensor(nCells=10)
  list[nDrug, nGene, nCell] = dim(T1)
 
  # Add increasing noise per drug
  T2 = T1
  for(d in 1:nDrug){
    T2[d,,] = T2[d,,] + rnorm(n=nGene*nCell, mean=0, sd=d*(1/nDrug))
  }
  out = ComputePCT_AllModes(T1, list(tmp = T2))
  y = smooth.spline(out$PCTd$tmp)$y
  stopifnot(identical(sort(y, decreasing=TRUE), y))
  
  # Add increasing noise per gene
  T2 = T1
  for(g in 1:nGene){
    T2[,g,] = T2[,g,] + rnorm(n=nDrug*nCell, mean=0, sd=g*(1/nGene))
  }
  out = ComputePCT_AllModes(T1, list(tmp = T2))
  y = smooth.spline(out$PCTg$tmp)$y
  stopifnot(identical(sort(y, decreasing=TRUE), y))
}

TestComputePCTPerSig = function(){
  # Make two random tensors, and set one equal to another in two locations
  set.seed(123)
  T1 = MakeTensor()
  T2 = MakeTensor()

  T1[5,,1] = T2[5,,1]
  T1[10,,2] = T2[10,,2]
  T1[4,,1] = T2[10,,2]
  
  # Make sure I can recover those locations
  PCTs = ComputePCTPerSig(T1, T2, format='list')
  
  stopifnot(PCTs$R[5,1] == 1)
  stopifnot(PCTs$R[10,2] == 1)
  stopifnot(PCTs$R[4,1] < 1)
}

TestTensor2Vec = function(){
  T1 = MakeTensor(NA_frac=0.1)
  T2 = MakeTensor(NA_frac=0.1)
  nEntries = prod(dim(T1))
  list[x1,x2] = Tensor2Vec(T1, T2)
  stopifnot(abs(length(x1) - nEntries*0.81) < (nEntries/100) )
}

TestComputeCSpPres = function(){}

#### Helper functions ####

MakeTensor = function(nDrugs=500, nGenes=10, nCells=2, NA_frac=0, fill=NA, noise=0, removeSig_frac=0){
  nEntries = nDrugs * nGenes * nCells
  
  nm = list(drugs=paste0('drug.', AlphaNames(nDrugs)),
            genes=paste0('gene.', AlphaNames(nGenes)),
            cells=paste0('cell.', AlphaNames(nCells)))

  if(is.na(fill)){
    tensor = array(data=rnorm(nEntries), dim=c(nDrugs,nGenes,nCells), dimnames=nm)
  }else{
    tensor = array(data=fill, dim=c(nDrugs,nGenes,nCells), dimnames=nm)
  }

  if(noise > 0){
    tensor = tensor + array(data=rnorm(nEntries, mean=0, sd=noise), dim=c(nDrugs, nGenes, nCells), dimnames=nm)
  }

  if(NA_frac > 0 && (removeSig_frac > 0|!is.na(fill))){
    stop('NA_frac is not compatible with other arguments')
  }

  if(NA_frac > 0){
    stopifnot(NA_frac <= 1)
    tensor = RemoveRandomEltsFromTensor(tensor, NA_frac=NA_frac)
  }else if(removeSig_frac > 0){
    stopifnot(removeSig_frac <= 1)
    idx = sample(nDrugs*nCells, size=nDrugs*nCells*removeSig_frac, replace=FALSE)
    toRemove = melt(tensor[,1,])[idx,]
    for(i in 1:nrow(toRemove)){
      tensor[toRemove$drugs[i],,toRemove$cells[i]] = NA
    }
  }
  
  return(tensor)
}

RemoveRandomEltsFromTensor = function(tensor, NA_frac=0.1){
  a = sample(tensor, size=prod(dim(tensor))*NA_frac)
  tensor[tensor %in% a] = NA
  return(tensor)
}

