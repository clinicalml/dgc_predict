TestSubsetTensorBy = function(){
  warning('add test')
}

TestGetDrugSlice = function(){
  warning('add test')
}

TestGetLincsAnnot = function(){
  L = GetLincsAnnot()
  stopifnot(nrow(L) == 20452)
  stopifnot(ncol(L) == 28)
}

TestNormSigs = function(){
  T1 = array(data=rnorm(120), dim=c(3,10,4))
  T2 = NormSigs(T1)
  T3 = NormSigs(T2)
  stopifnot(max(abs(T2-T3)) < 1e-8)
  stopifnot(Norm2(T3[1,,1]) == 1)
}

TestComputeDensity = function(){
  load(DataDir('tensors/T_test.RData'))
  d1 = length(which(!is.na(T_test))) / prod(dim(T_test))
  d2 = ComputeDensity(T_test)
  stopifnot(d1 == d2)
}

TestGetLmGenes = function(){
  nGene = 978 

  lmGenes = GetLmGenes('symbol')
  stopifnot(length(lmGenes) == nGene)
  stopifnot(class(lmGenes) == 'character')
  stopifnot(lmGenes[1] == 'AARS')
  
  lmGenes = GetLmGenes('entrez')
  stopifnot(length(lmGenes) == nGene)
  stopifnot(class(lmGenes) == 'integer')
  stopifnot(lmGenes[1] == 16)
  
  lmGenes = GetLmGenes('probe')
  stopifnot(length(lmGenes) == nGene)
  stopifnot(class(lmGenes) == 'character')
  stopifnot(lmGenes[1] == '201000_at')
}

TestMapEntrez2Uniprot = function(){
  x = MapEntrez2Uniprot()
  stopifnot(names(x)[1] == '1')
  stopifnot(x[[1]][1] == 'P04217')
}

TestMapEntrez2Hugo = function(){
  stopifnot(MapEntrez2Hugo(55603) == 'FAM46A')
  stopifnot(MapEntrez2Hugo('55603') == 'FAM46A')
  stopifnot(all(MapEntrez2Hugo(c(100302217,23559)) == c('MIR1827', 'WBP1')))
  stopifnot(all(MapEntrez2Hugo(list(100302217,23559)) == c('MIR1827', 'WBP1')))
}

TestMapUniprot2Entrez = function(){
  proteins = c('F5GWI4','A0AAB7', 'B6EC88')
  entrez = MapUniprot2Entrez(proteins, printFlag=FALSE)
  
  stopifnot(length(entrez) == length(proteins))
  mapBack = MapEntrez2Uniprot()
  stopifnot(proteins[c(1,3)] %in% unlist(mapBack[as.character(entrez)]))
}

TestGetLincs2Pubchem = function(){
  map = GetLincs2Pubchem()
  stopifnot(dim(map) == c(20452,2))
  stopifnot(class(map$pert_id) == 'character')
  stopifnot(class(map$pubchem_id) == 'character')
  stopifnot(length(unique(map$pert_id)) > 1000)
  stopifnot(length(unique(map$pubchem_id)) > 1000)
}

TestDataFile = function(){
  file = DataFile('tensors/cmap.mat')
  stopifnot(file.exists(file))
}

TestNumSigs = function(){
  load(DataDir('tensors/T_small.RData'))
  n = NumSigs(T_small)
  stopifnot(n == 36)
  
  nd = NumSigs(T_small, dim='drug')
  nc = NumSigs(T_small, dim='cell')
  stopifnot(sum(nd) == n)
  stopifnot(sum(nc) == n)
  stopifnot(length(nd) == dim(T_small)[1])
  stopifnot(length(nc) == dim(T_small)[3])
}

TestGetTensorAnnot = function(){
  annot = GetTensorAnnot()
  stopifnot(c('geneIds','cellIds','pertIds', 
              'pubchemIds', 'geneSymbols', 'tissues',
              'targets_entrez', 'targets', 'pertName',
              'pertSummary') %in% names(annot))
  
  d = GetTensorDims()
  stopifnot(length(annot$geneIds) == d$nGene)
  stopifnot(length(annot$cellIds) == d$nCell)
  stopifnot(length(annot$pertIds) == d$nDrug)
  stopifnot(length(annot$pubchemIds) == d$nDrug)
  stopifnot(length(annot$geneSymbols) == d$nGene)
  stopifnot(length(annot$tissues) == d$nCell)
  stopifnot(length(annot$targets_entrez) == d$nDrug)
  stopifnot(length(annot$targets) == d$nDrug)
  stopifnot(length(annot$pertName) == d$nDrug)
  stopifnot(length(annot$pertSummary) == d$nDrug)
  
  for(a in names(annot)){stopifnot(class(annot[[a]]) == 'character')}

  stopifnot(annot$pertName[599] == 'triciribine')
  stopifnot(annot$targets[599] == 'AKT1')
  stopifnot(annot$geneIds[10] == '6657')
  stopifnot(annot$geneSymbols[10] == 'SOX2')
}

TestUnfoldTensor = function(){
  load(DataDir('tensors/T_small.RData'))
  X = T_small
  M1 = UnfoldTensor(X, 1)
  M2 = UnfoldTensor(X, 2)
  M3 = UnfoldTensor(X, 3)
  stopifnot(rownames(M1) == dimnames(X)[[1]])
  stopifnot(rownames(M2) == dimnames(X)[[2]])
  stopifnot(rownames(M3) == dimnames(X)[[3]])
  
  n = dim(X)[1]
  list[a, b] = unlist(strsplit(colnames(M1)[n], '[.]'))
  stopifnot(X[1,a,b] == M1[1,n])
  
  list[a, b] = unlist(strsplit(colnames(M2)[n], '[.]'))
  stopifnot(X[a,1,b] == M2[1,n])
  
  list[a, b] = unlist(strsplit(colnames(M3)[n], '[.]'))
  stopifnot(X[a,b,1] == M3[1,n])
}

TestTensorZ = function(){
  load(DataDir('tensors/T_small.RData'))
  X = T_small
  Z = TensorZ(X)$Z
  m = apply(Z, 2, function(x) mean(x, na.rm=TRUE))
  v = apply(Z, 2, function(x) var(as.vector(x), na.rm=TRUE))
  stopifnot(Norm2(m) < 1e-10)
  stopifnot(Norm2(v - 1) < 1e-10)
  stopifnot(max(abs((Z[,,1] - TensorZ(Z)$Z[,,1])), na.rm=TRUE) < 1e-10)
  
  for(i in 1:4){
    x = as.vector(X[,i,])
    z = as.vector(Z[,i,])
    stopifnot(cor(x,z,use='pairwise.complete')==1)
  }
}

TestTensorDEG = function(){
  load(DataDir('tensors/T50_1.RData'))
  X = T_meas
  eps = 1
    for(perc in c(10,20)){
      for(sym in c(TRUE, FALSE)){
        out = TensorDEG(X, normGene=FALSE, method='sig', percDEG=perc, symmetric=sym)
        stopifnot(abs(perc-out$p) < eps)
        CheckXD(x = X[2,,2], d = out$D[2,,2])
      }
    }
}
  
CheckXD = function(x,d){
  x = na.omit(as.vector(x))
  d = na.omit(as.vector(d))
  stopifnot(length(x)>0)
  stopifnot(d[which(x==min(x))] == -1)
  stopifnot(d[which(x==max(x))] == 1)
  stopifnot(d[which(x == min(abs(x)))] == 0)
}

TestCallDEG = function(){
  n = 10000
  x = array(data=rnorm(n), dim=c(100,100))
  x[sample(1:n,100)] = NA
  nn = length(which(!is.na(x)))
  p = 0.99
  
  for(sym in c(TRUE, FALSE)){
    d = CallDEG(x, p, symmetric=sym)
    stopifnot(identical(sort(unique(as.vector(d))), c(-1,0,1)))
    noDEGs = length(which(abs(d) == 0))
    stopifnot(abs(p - (noDEGs/nn)) <= 0.005)
    top10 = order(x, decreasing=TRUE)[1:10]
    bottom10 = order(x)[1:10]
    stopifnot(all(d[top10] == 1))
    stopifnot(all(d[bottom10] == -1))
    nDEGs = length(which(abs(d)==1))
    expDEGs = round(n*(1-p))
    stopifnot(abs(nDEGs - expDEGs) <= 2)
  }
}

TestComputeAUC = function(){
  n = 1000
  est = c(rep(5,n), rep(5.1,n))
  labels = c(rep(0,n), rep(1, n))
  stopifnot(ComputeAUC(est, labels)==1)
  
  labels = sample(labels)
  stopifnot(abs(ComputeAUC(est, labels)-0.5) < 0.05)
}

TestComputeCellSpecificity = function(){
  # construct one tensor with identical values
  T1 = array(data=1, dim = c(2,20,5))
  cs1 = ComputeCellSpecificity(T1, normalize=TRUE)$cs
  stopifnot(all(abs(cs1) < 1e-10))
  
  # construct another tensor with random values
  T2 = array(data=rnorm(1000), dim=c(2,100,5))
  cs2 = ComputeCellSpecificity(T2, normalize=TRUE)$cs
  stopifnot(all(abs(cs2-1) < 0.1))
}

TestGetCellSpecificity = function(){
  cs = GetCellSpecificity()
  stopifnot(all(names(cs) == GetTensorAnnot()$pertIds))
  stopifnot(all(cs > 0 && cs < 2))
}

TestComputeGeneGeneCor = function(){
  load(DataDir('tensors/T50_1.RData'))
  X = T_meas
  G = ComputeGeneGeneCor(X, nGene=10, cellSpecific=FALSE)
  gList1 = list()
  for(i in 1:dim(X)[3]){
    gList1[[i]] = G
  }
  gList2 = ComputeGeneGeneCor(X, nGene=10, cellSpecific=TRUE, print=FALSE)
  cm = CorMatrixList(gList1, gList2)
  stopifnot(all(cm > 0.2))
}

TestGetGeneGeneCor = function(){
  G = GetGeneGeneCor()
  nGene = GetTensorDims('gene')
  stopifnot(dim(G) == c(nGene, nGene))
  stopifnot(IsSymmetric(G))
  stopifnot(all(diag(G) == 1))
  stopifnot(min(G) >= -1)
  stopifnot(max(G) == 1)
  stopifnot(any(is.na(G)) == FALSE) 
}

TestRankSigs = function(){
  load(DataDir('tensors/T50_1.RData'))
  X = abs(T_meas[1:20,1:50,])
  cs = ComputeCellSpecificity(X)$cs
  R = RankSigs(X)
  csR = ComputeCellSpecificity(R)$cs
  stopifnot(range(R, na.rm=TRUE) == c(1, 50))
  stopifnot(cor(cs, csR) > 0.4)
  stopifnot(which.max(X[1,,5]) == which.max(R[1,,5]))
}

TestGetAccuracyPerMode = function(){
  C = GetAccuracyPerMode()
  methods = c('mean', 'mean2', 'knn', 'tensor')
  modes = c('gene', 'drug', 'cell')
  stopifnot(names(C) == methods)
  stopifnot(all(modes %in% names(C$mean)))
  for(method in methods){
    for(mode in modes){
      r = range(C[[method]][[mode]])
      stopifnot(r[1] > -0.1)
      stopifnot(r[2] < 1)
    }
  }
  stopifnot(cor(C$mean$gene, C$mean2$gene) > 0.5)
  stopifnot(cor(C$tensor$gene, C$mean2$gene) > 0.5)
  stopifnot(cor(C$tensor$drug, C$mean2$drug) > 0.5)
}

TestGetTensorDims = function(){
  a = GetTensorDims()
  T = GetDataTensor()
  stopifnot(dim(T)[1] == a$nDrug)
  stopifnot(dim(T)[2] == a$nGene)
  stopifnot(dim(T)[3] == a$nCell)
  
  nDrug = GetTensorDims('drug')
  nGene = GetTensorDims('gene')
  nCell = GetTensorDims('cell')
  
  stopifnot(a$nDrug == nDrug)
  stopifnot(a$nGene == nGene)
  stopifnot(a$nCell == nCell)
}

TestGetExpCount = function(){
  expCount = GetExpCount()
  annot = GetTensorAnnot()
  stopifnot(all(rownames(expCount) == annot$pertIds))
  stopifnot(all(colnames(expCount) == annot$cellIds))
  stopifnot(min(expCount) == 0)
  stopifnot(all(!is.na(expCount)))
}

TestCheckColumnStructure = function(){
  load(DataDir('tensors/T_test.RData'))
  stopifnot(CheckColumnStructure(T_test))
  
  T2 = T_test
  T2[1,2,1] = NA
  stopifnot(!CheckColumnStructure(T2))
  
  T3 = T_test
  T3[1,1,1] = NA
  stopifnot(!CheckColumnStructure(T3))
  
  T3[5,1,2] = 5
  stopifnot(!CheckColumnStructure(T3))
}

TestRemoveMissingDataFromXY = function(){
  X = matrix(data=1:12, nrow=6, ncol=2)
  Y = matrix(data=1:6, nrow=6, ncol=1)
  
  list[X1,Y1] = RemoveMissingDataFromXY(X,Y)
  stopifnot(identical(X,X1))
  stopifnot(identical(Y,Y1))
  
  X[2,] = NA
  X[3,1] = NA
  
  list[X1,Y1] = RemoveMissingDataFromXY(X,Y)
  stopifnot(identical(X[-2:-3,],X1))
  stopifnot(identical(Y[-2:-3,],Y1))
  
  Y[4] = NA
  
  list[X1,Y1] = RemoveMissingDataFromXY(X,Y)
  stopifnot(identical(X[-2:-4,],X1))
  stopifnot(identical(Y[-2:-4,],Y1))
}

TestLoadTensorMat = function(){}

TestGetGeneIdsTensor = function(){
  geneIds = GetGeneIdsTensor()
  stopifnot(length(geneIds)==978)
  stopifnot(is.character(geneIds))
  stopifnot(geneIds[1] == '5720')
  stopifnot(geneIds[978] == '874')
}
