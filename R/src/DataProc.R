
SubsetTensorDims = function(tensorWithDesiredDims, tensorToSubset){
  return(tensorToSubset[dimnames(tensorWithDesiredDims)[[1]], 
                        dimnames(tensorWithDesiredDims)[[2]],
                        dimnames(tensorWithDesiredDims)[[3]]])
}

GetDrugSlice = function(tensor, drug){
  return(t(na.omit(t(tensor[drug,,]))))
}

GetLincsAnnot = function(){
  load(DataDir('metadata/lincsAnnot.RData'))
  return(L)
}

NormSigs = function(A){
  for(i in 1:nrow(A)){
    for(j in 1:dim(A)[3]){
      A[i,,j] = A[i,,j] / Norm2(A[i,,j])
    }
  }
  return(A)
}

ComputeDensity = function(tensor){
  numPresent = NumSigs(tensor)
  numPossible = dim(tensor)[1]*dim(tensor)[3]
  return(numPresent/numPossible)
}

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

MapEntrez2Uniprot = function(){
  x = org.Hs.egUNIPROT
  mapped_genes = mappedkeys(x) 
  return(as.list(x[mapped_genes]))
}

MapEntrez2Hugo = function(entrez_ids){
  if(!all(is.na(entrez_ids))){
    load(DataDir('metadata/hgnc_to_entrez.RData'))
    names(gene_annot) = c('hugo_id', 'entrez_id')
    idx_remove = which(gene_annot$hugo_id == '')
    gene_annot = gene_annot[-idx_remove,]
    idx = match(entrez_ids, gene_annot$entrez_id)
    out = gene_annot$hugo_id[idx]
  }else{
    out = rep(NA, length(entrez_ids))
  }
  return(out)
}

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

GetLincs2Pubchem = function(){
  load(DataFile('metadata/lincsAnnot_inHouse_pertInfo_allCompounds.RData'))
  lincsIdMap = lincsAnnot[,c('pert_id', 'pubchem_cid.pertInfo')]
  names(lincsIdMap) = c('pert_id', 'pubchem_id')
  return(lincsIdMap)
}

DataFile = function(file){
  return(paste0(GetConfig('DATAPATH'),'/',file))
}

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

UnfoldTensor = function(X, dim=1){
  return(wrap(X, map=list(dim,NA)))
}

# standardize tensor values per gene
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

# x is a gene expression profile
# p is the percentile threshold (between 0 and 1)
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

ComputeAUC = function(est, labels, computeROC=FALSE, abs=TRUE){
  if(abs){est = abs(est)}
  pred = prediction(est,labels)
  auc = performance(pred, measure = "auc")
  if(computeROC){
    perf = performance(pred, measure = "tpr", x.measure = "fpr")
    roc = data.frame(fpr=unlist(perf@x.values),tpr=unlist(perf@y.values))
    out = list(auc=auc@y.values[[1]], roc=roc)
  }else{
    out = auc@y.values[[1]]
  }
  return(out)
}

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

LoadTensorMat = function(file){
  out = readMat(file)
  tensor = out$T
  geneIds = unlist(out$geneIds)
  pertIds = unlist(out$pertIds)
  cellIds = unlist(out$cellIds)
  dimnames(tensor) = list(pertIds, geneIds, cellIds)
  return(list(tensor=tensor, pertIds=pertIds, geneIds=geneIds, cellIds=cellIds))
}

GetGeneIdsTensor = function(){
  load(DataDir('tensors/T50_1.RData'))
  return(dimnames(T_meas)[[2]])
}
