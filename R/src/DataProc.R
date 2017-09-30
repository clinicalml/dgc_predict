library(ROCR)

# Restrict/reorder the dimensions of the second argument to those of the first
# argument. Both should be 3-dimensional tensors, and it is assumed that the
# second tensor contains at least all of the dimensions of the first.
SubsetTensorDims = function(tensorWithDesiredDims, tensorToSubset){
  return(tensorToSubset[dimnames(tensorWithDesiredDims)[[1]], 
                        dimnames(tensorWithDesiredDims)[[2]],
                        dimnames(tensorWithDesiredDims)[[3]]])
}

# Extract one slice of the tensor corresponding to a chosen drug. The drug
# reference can be named or just an index.
GetDrugSlice = function(tensor, drug){
  return(t(na.omit(t(tensor[drug,,]))))
}

# Annotation and metadata (e.g. targets, side effects) for all compounds in LINCS L1000
GetLincsAnnot = function(){
  load(DataDir('metadata/lincsAnnot.RData'))
  return(L)
}

# Normalize each gene profile (or the second dimension) to have L2 norm of 1
NormSigs = function(A){
  for(i in 1:nrow(A)){
    for(j in 1:dim(A)[3]){
      A[i,,j] = A[i,,j] / Norm2(A[i,,j])
    }
  }
  return(A)
}

# Observation density in the tensor. Assumes the pattern of missing entries is
# consistent across all genes, i.e. each drug/cell combination is either
# completely observed or completely absent.
ComputeDensity = function(tensor){
  numPresent = NumSigs(tensor)
  numPossible = dim(tensor)[1]*dim(tensor)[3]
  return(numPresent/numPossible)
}

# Gene information for 978 landmark genes
GetLmGenes = function(type='symbol'){
  file = DataFile('metadata/L1000_CD_landmark_geneInfo.txt')
  lmGeneInfo = Factor2Char(read.table(file, header=TRUE, sep='\t'))
  if(type == 'symbol'){
    out = lmGeneInfo$gene_symbol
  }else if(type == 'entrez'){
    out = lmGeneInfo$gene_id
  }else if(type == 'probe'){
    out = lmGeneInfo$probe_id
  }else if(type == 'all'){
    out = lmGeneInfo
  }else{
    warning('unexpected type, returning entrez id')
    out = lmGeneInfo$gene_id
  }
  return(out)
}

# Map gene entrez ID to Uniprot ID
MapEntrez2Uniprot = function(){
  x = org.Hs.egUNIPROT
  mapped_genes = mappedkeys(x) 
  return(as.list(x[mapped_genes]))
}

# Map gene entrez ID to Hugo ID
MapEntrez2Hugo = function(entrez_ids){
  if(!all(is.na(entrez_ids))){
    load(DataDir('metadata/hgnc_to_entrez.RData'))
    names(map) = c('hugo_id', 'entrez_id')
    idx_remove = which(map$hugo_id == '')
    map = map[-idx_remove,]
    idx = match(entrez_ids, map$entrez_id)
    out = map$hugo_id[idx]
  }else{
    out = rep(NA, length(entrez_ids))
  }
  return(out)
}

# Map Entrez ID to HGNC symbol
MapEntrezToSymbol = function(entrez_ids, lm, map=NA){
  if(is.na(map)){
    if(lm){
      map = GetLmGenes('all')
    }else{
      load(DataDir('metadata/hgnc_to_entrez.RData'))
    }
  }
  stopifnot(c('gene_id', 'gene_symbol') %in% names(map))
  idx = match(entrez_ids, map$gene_id)
  if(all(is.na(idx))){
    warning('no ids matched')
  }
  return(map$gene_symbol[idx])
}

# Map Uniprot id to entrez
MapUniprot2Entrez = function(proteins, printFlag=TRUE, collapse=FALSE){
  map = MapEntrez2Uniprot()
  p2e = list()
  
  for(p in proteins){
    p2e[[p]] = NA
    if(printFlag) print(p)
    for(i in 1:length(map)){
      if(p %in% map[[i]]){
        p2e[[p]] = c(p2e[[p]], names(map)[[i]]) 
      }
    }
    if(length(p2e[[p]]) >= 2){
      p2e[[p]] = setdiff(p2e[[p]], NA)
      if(collapse){
        p2e[[p]] = paste(p2e[[p]], collapse='|')
      }
    }
  }
  return(as.character(p2e))
}

# Get full path to file under data directory
DataFile = function(file){
  return(paste0(GetConfig('DATAPATH'),'/',file))
}

# Compute the number of observed profiels in the tensor, assuming the
# column-structured pattern of missing entries
NumSigs = function(T, dim='all'){
  if(dim == 'all'){
    A = T[,1,]
    out = length(which(!is.na(A)))
  }else if(dim == 'drug'){
    out = apply(T[,1,], 1, function(x) length(which(!is.na(x))))
  }else if(dim == 'cell'){
    out = apply(T[,1,], 2, function(x) length(which(!is.na(x))))
  }else{
    error('unexpected argument for dim')
  }
  return(out)
}

# Unfold the tensor along the chosen dimension to generate a matrix.
UnfoldTensor = function(X, dim=1){
  return(wrap(X, map=list(dim,NA)))
}

# Standardize tensor values per gene
TensorZ = function(X, m=NULL, v=NULL){
  if(is.null(m)){
    m = apply(X, 2, function(x) mean(x, na.rm=TRUE))
  }
  if(is.null(v)){
    v = apply(X, 2, function(x) var(as.vector(x), na.rm=TRUE))
  }
  Z = array(data=NA, dim=dim(X), dimnames=dimnames(X))
  for(i in 1:length(m)){
    Z[,i,] = (X[,i,] - m[i]) / sqrt(v[i])
  }
  return(list(Z=Z,m=m,v=v))
}

# Takes input tensor X and labels each observation as differentially expressed
# or not (1 = up-reg, -1=down, 0=unchanged). Output p is the percent of
# observations labeled as DEG
TensorDEG = function(X, normGene=FALSE, method='tensor', percDEG=2, symmetric=FALSE){
  if(all(c('drug', 'gene', 'cell') %in% names(dimnames(X)))){
    X = aperm(X, c('drug','gene', 'cell'))
  }

  if(normGene){# normalize each gene's values across tensor
    Z = TensorZ(X)$Z
  }else{
    Z = X
  }
  
  # define a few variables
  p = 1- percDEG/100

  # initilize output
  D = array(data=0, dim=dim(Z), dimnames=dimnames(Z))
  D[is.na(Z)] = NA
  
  # identify DEGs
  if(method == 'tensor'){
    D = CallDEG(Z, p, symmetric=symmetric)
  }else if(method == 'sig'){
    for(drug in 1:dim(D)[1]){
      for(cell in 1:dim(D)[3]){
        D[drug,,cell] = CallDEG(Z[drug,,cell], p, symmetric=symmetric)
      }
    }
  }else if(method == 'gene'){
    for(gene in 1:dim(D)[2]){
      D[,gene,] = CallDEG(Z[,gene,], p, symmetric=symmetric)
    }
  }else{
    stop('no DEGs called, unexpected method')
  }
  p_actual = 100*length(which(abs(D)==1)) / length(which(!is.na(D)))
  return(list(D=D, p=p_actual))
}

# x is a gene expression profile (i.e. a vector) p is the percentile threshold
# (between 0 and 1). Larger p corresponds to fewer genes being labeled as DEGs.
CallDEG = function(x, p, symmetric=FALSE){
  d = x
  d[which(!is.na(x))] = 0
  n_elts = length(which(!is.na(x)))
  if(length(which(d == 0) > 0)){
    if(symmetric){
      n = round( (1-p)*n_elts / 2 )
      if(n > 0){
        d[order(x)[1:n]] = -1
        d[order(x, decreasing=TRUE)[1:n]] = 1
      }
    }else{
      q = quantile(abs(x), probs=p, na.rm=TRUE)
      d[which(x >= q)] = 1
      d[which(x <= -q)] = -1
    } 
  }
  return(d)
}

# Compute AUC from some estimated scores in 'est', comparing to the true labels
# (should be the same length), Can also compute ROC curves to be plotted.
ComputeAUC = function(est, labels, computeROC=FALSE, abs=FALSE, na.rm=FALSE){
  if(na.rm){
    idxNA = which(is.na(est))
    idxLab = which(is.na(labels))
    idxRemove = union(idxNA, idxLab)
    if(length(idxRemove) > 0){
      print(sprintf('  removing %d missing values in ComputeAUC', length(idxRemove)))
      est = est[-idxRemove]
      labels = labels[-idxRemove]
    }
  }
  if(length(unique(labels)) != 2){
    warning(sprintf('labels have %d unique entries, should have 2', length(unique(labels))))
    out = NA
  }else{
    if(abs){
      est = abs(est)
    }
    pred = prediction(est,labels)
    auc = performance(pred, measure = "auc")
    if(computeROC){
      perf = performance(pred, measure = "tpr", x.measure = "fpr")
      roc = data.frame(fpr=unlist(perf@x.values),tpr=unlist(perf@y.values))
      out = list(auc=auc@y.values[[1]], roc=roc)
    }else{
      out = auc@y.values[[1]]
    }
  }
  return(out)
}

# Defines cell specificity for each drug in tensor X, based on mean pairwise
# cosine distances, as described in the paper. If the signatures in the tensor
# are already normalized, you can save compute time by setting normalize to
# FALSE.
ComputeCellSpecificity = function(X, normalize=FALSE){
  nDrug = dim(X)[1]
  nCell = dim(X)[3]
  # compute pairwise cosine distances between drug signatures in measured tensor
  A = array(data=NA, dim=c(nCell,nCell))
  cs = rep(NA, times=nDrug)
  r = list()
  for(d in 1:nDrug){
    C = A
    for(c1 in 1:nCell){ 
      x = X[d,,c1]
      if(any(!is.na(x))){
        for(c2 in IncreasingSequence(c1+1,nCell)){
          y = X[d,,c2]
          if(any(!is.na(y))){
            C[c1,c2] = CosineDistance(x, y, normalize=normalize)
          }
        }
      }
    }
    cs[d] = mean(as.numeric(C), na.rm=TRUE)
    r[[d]] = range(as.numeric(C), na.rm=TRUE)
  }
  names(cs) = dimnames(X)[[1]] 
  names(r) = dimnames(X)[[1]]
  return(list(cs=cs, r=r))
}

# Compute gene gene correlation from tensor, either per cell type (if
# cellSpecific = TRUE), or across all cell types
ComputeGeneGeneCor = function(tensor, nGene=dim(tensor)[2], cellSpecific=FALSE, print=TRUE){
  geneIds = dimnames(tensor)[[2]][1:nGene]
  
  if(!cellSpecific){
    G = matrix(data=NA, nrow=nGene, ncol=nGene, dimnames=list(geneIds, geneIds))
    for(i in 1:nGene){
      g_i = as.vector(tensor[,i,])
      for(j in IncreasingSequence(i+1,nGene)){
        g_j = as.vector(tensor[,j,])
        G[i,j] = cor(g_i, g_j, use='pairwise')
      }
    }
    diag(G) = 1
    out = SymmetrifyMatrix(G, diag=NA)
  }else{
    cellIds = dimnames(tensor)[[3]]
    out = list()
    for(c in 1:length(cellIds)){
      cell = cellIds[c]
      if(print){print(cell)}
      out[[cell]] = cor(tensor[,1:nGene,cell], use='pairwise')
    }
  }
  return(out)
}

# Converts tensor values into ranked gene lists per gene expression profile
RankSigs = function(tensor){
  if(any(tensor < 0, na.rm=TRUE)){
    warning('Did you want to take absolute value?')
  }
  dims = dim(tensor)
  for(i in 1:dims[1]){
    for(j in 1:dims[3]){
      if(!is.na(tensor[i,1,j])){
        tensor[i,,j] = rank(tensor[i,,j]) # this gives the largest values the highest (i.e. largest valued) ranks
      }
    }
  }
  return(tensor)
}

# Verify that pattern of missing entries in the tensor is column-structured
CheckColumnStructure = function(tensor){
  columnStructured = TRUE
  A = is.na(tensor[,1,])
  for(i in 2:dim(tensor)[2]){
    B = is.na(tensor[,i,])
    if(!identical(A,B)){
      columnStructured = FALSE
      break;
    }
  }
  return(columnStructured)
}

# Load tensor from .mat file
LoadTensorMat = function(file){
  out = readMat(file)
  tensor = out$T
  geneIds = unlist(out$geneIds)
  pertIds = unlist(out$pertIds)
  cellIds = unlist(out$cellIds)
  dimnames(tensor) = list(pertIds, geneIds, cellIds)
  return(list(tensor=tensor, pertIds=pertIds, geneIds=geneIds, cellIds=cellIds))
}

# Load measured tensor as well as cross-validated tensor values for all four methods presented in the paper
LoadTensors = function(tsize='small', print=FALSE){
  library(rhdf5)
  tensors = list()
  
  if(print){print('Loading data tensor..')}
  
  tensors$meas = LoadTensorMat(DataDir(sprintf('tensors/%s.mat', tsize)))$tensor
  names(dimnames(tensors$meas)) = c('drug','gene','cell')
  
  if(print){print('Loading cross-validated tensors..')}
  
  file = ResultsDir(sprintf('%s/%s_tensor_cv_results.mat', tsize, tsize))
  
  tensors$cv = list(mean=h5read(file, '#refs#/b'), mean2=h5read(file,'#refs#/c'),
                    dnpp=h5read(file,'#refs#/d'), tensor=h5read(file,'#refs#/e'))
  
  tensors$cv = lapply(tensors$cv, function(tensor){
    dimnames(tensor) = dimnames(tensors$meas); return(tensor)})
  
  return(tensors)
}

# Get the landmark Gene ID ordering in the tensor (consistent for all tensors)
GetGeneIdsTensor = function(){
  load(DataDir('tensors/test/T50_1.RData'))
  return(dimnames(T_meas)[[2]])
}

# Get vector of method names
CompletionMethods = function(version='matlab', knn=FALSE){
  if(version == 'matlab'){
    if(knn){methods=c('mean','mean2', 'knnd', 'fa_lrtc')}
    else{methods=c('mean','mean2', 'dnpp', 'fa_lrtc')}
  }else if(version == 'R'){
    if(knn){methods=c('mean','mean2', 'knnd', 'tensor')}
    else{methods=c('mean','mean2', 'dnpp', 'tensor')}
  }else if(version == 'title'){
    if(knn){methods = c('1D-Mean', '2D-Mean', 'KNN', 'Tensor')}
    else{methods = c('1D-Mean', '2D-Mean', 'DNPP', 'Tensor')}
  }else{
    stop('Unexpected value for version')
  }
  return(methods)
}
