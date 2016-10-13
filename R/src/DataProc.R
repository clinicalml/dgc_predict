library(org.Hs.eg.db)
library(R.utils)
library(plotrix)
library(ROCR)

GetLincsAnnot <- function(){
  load(DataDir('annot/lincs_annot/pertInfo/lincsAnnot_v2.RData'))
  return(L)
}

NormSigs <- function(A){
  for(i in 1:nrow(A)){
    for(j in 1:dim(A)[3]){
      A[i,,j] = A[i,,j] / Norm2(A[i,,j])
    }
  }
  return(A)
}

ComputeDensity <- function(tensor){
  numPres = NumSigs(tensor)
  numPos = dim(tensor)[1]*dim(tensor)[3]
  return(numPres/numPos)
}

GetLmGenes <- function(type='symbol'){
  file <- DataFile('annot/lincs_annot/geneInfo/L1000_CD_landmark_geneInfo.txt')
  lmGeneInfo <- Factor2Char(read.table(file, header=TRUE, sep='\t'))

  if(type == 'symbol'){
    out <- lmGeneInfo$gene_symbol
  }else if(type == 'entrez'){
    out <- lmGeneInfo$gene_id
  }else if(type == 'probe'){
    out <- lmGeneInfo$probe_id
  }else if(type == 'all'){
    out <- lmGeneInfo
  }else{
    warning('unexpected type, returning entrez id')
    out <- lmGeneInfo$gene_id
  }
  return(out)
}

CombineBinaryMatrices <- function(M1, M2){
  M1 <- as.matrix(M1)
  M2 <- as.matrix(M2)
  
  r1 <- rownames(M1)
  r2 <- rownames(M2)
  c1 <- colnames(M1)
  c2 <- colnames(M2)
  
  stopifnot(length(intersect(r1, r2)) > 0)
  stopifnot(length(intersect(c1, c2)) > 0)
  
  allRows <- union(rownames(M1), rownames(M2))
  allCols <- union(colnames(M1), colnames(M2))
 
  M <- matrix(data=0, nrow=length(allRows), ncol=length(allCols), 
              dimnames=list(allRows, allCols))
  M[r1,c1] <- M1
  M[r2,c2] <- MatrixCast(M[r2,c2] | M2, 'numeric')

  return(M)
}

# input should be two column data frame representing an edge list
Df2BinaryMatrix <- function(df){
  stopifnot(ncol(df) == 2)
  R <- unique(as.character(df[,1]))
  C <- unique(as.character(df[,2]))
  M <- matrix(data=0, nrow=length(R), ncol=length(C), dimnames=list(R,C))
  for(i in 1:nrow(df)){
    M[df[i,1],df[i,2]] <- 1
  }
  return(M)
}

# input should be a 2 or 3 column data frame, representing an edge list and optionally edge weights
Df2Matrix <- function(df){
  if(ncol(df) == 2){
    M <- Df2BinaryMatrix(df)
  }else if(ncol(df) == 3){
    R <- unique(as.character(df[,1]))
    C <- unique(as.character(df[,2]))
    M <- matrix(data=0, nrow=length(R), ncol=length(C), dimnames=list(R,C))
    for(i in 1:nrow(df)){
      M[df[i,1],df[i,2]] <- df[i,3]
    }
  }else{
    stop('unexpected number of columns')
  }
  return(M)
}

BinaryMatrix2Edges <- function(M){
  n <- sum(M)
  E <- Factor2Char(data.frame(from=rep('',n), to=rep('',n)))
  count <- 1
  for(from in rownames(M)){
    for(to in colnames(M)){
      if(M[from,to] == 1){
        E[count,] <- c(from, to)
        count <- count+1
      }
    }
  }
  return(E)
}

RestrictToLm <- function(genes){
  lmGenes <- GetLmGenes('entrez')
  for(i in 1:length(genes)){
    genes[[i]] <- genes[[i]][lmGenes %in% genes[[i]]]
    if(length(genes[[i]]) == 0){
      warning(sprintf('no genes in set %d in landmark set', i))
    }
  }
  return(genes)
}

MapEntrez2Uniprot <- function(){
  x <- org.Hs.egUNIPROT
  mapped_genes <- mappedkeys(x) 
  return(as.list(x[mapped_genes]))
}

MapEntrez2Hugo <- function(entrez_ids){
  if(!all(is.na(entrez_ids))){
    load(DataDir('/annot/gene_annot/hgnc_to_entrez.RData'))
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

MapUniprot2Entrez <- function(proteins, printFlag=TRUE, collapse=FALSE){
  map <- MapEntrez2Uniprot()
  p2e <- list()
  
  for(p in proteins){
    p2e[[p]] <- NA
    if(printFlag) print(p)
    for(i in 1:length(map)){
      if(p %in% map[[i]]){
        p2e[[p]] <- c(p2e[[p]], names(map)[[i]]) 
      }
    }
    if(length(p2e[[p]]) >= 2){
      p2e[[p]] <- setdiff(p2e[[p]], NA)
      if(collapse){
        p2e[[p]] <- paste(p2e[[p]], collapse='|')
      }
    }
  }
  return(as.character(p2e))
}

# GetCells2Keep <- function(){
#   return(c('A375','A549','ASC','HA1E','HCC515','HEPG2','HT29','MCF7','NPC','PC3','VCAP'))
# }

GetLincs2Pubchem <- function(){
  load(DataFile('/annot/lincs_annot/pertInfo/lincsAnnot_inHouse_pertInfo_allCompounds.RData'))
  lincsIdMap <- lincsAnnot[,c('pert_id', 'pubchem_cid.pertInfo')]
  names(lincsIdMap) <- c('pert_id', 'pubchem_id')
  return(lincsIdMap)
}

ComputeStringScore <- function(P){
  f <- 1000
  return(1-(1-(P$experimental)/f)*(1-(P$database)/f)*(1-(P$textmining)/f))
}

Str2GeneVec <- function(a){
  b <- gsub(pattern='"', replacement='', a)
  c <- unlist(strsplit(b, split=', '))
  if(length(c) > 0){
    c[1] <- substr(c[1], 3, nchar(c[1]))
    d <- c[length(c)]
    c[length(c)] <- substr(d, 1, nchar(d)-2)
  }
  return(c)
}

ConstructLmVec <- function(upGenes, downGenes, lmGenes=lmGenes){
  upBinary <- as.integer(lmGenes %in% upGenes)
  downBinary <- -(as.integer(lmGenes %in% downGenes))
  vec <- upBinary + downBinary
  if(sum(abs(vec))==0){
    browser()
  }
  return(upBinary + downBinary)
}

ExpandGraph <- function(currentV, currentE=NULL, remainingE, printFlag=FALSE){
  expandedV <- currentV
  
  if(is.null(currentE)){
    expandedE <- data.frame(V1=character(), V2=character())
  }else{
    expandedE <- currentE
  }
  
  for(v in currentV){
    idx1 <- which(remainingE$V1 == v)
    idx2 <- which(remainingE$V2 == v)
    idx <- union(idx1, idx2)
    
    if(length(idx) > 0){  
      # expand vertex and edge sets
      expandedV <- union(expandedV, c(remainingE$V2[idx1], remainingE$V1[idx2]))
      expandedE <- rbind(expandedE, remainingE[idx, ])
      
      # remove corresponding edges from remainder set
      remainingE <- remainingE[-idx,]
      
      if(printFlag){
        print(sprintf('%s added %d new edges', v, length(idx))) 
      }
    }
    
    if(nrow(remainingE) == 0){
      if(printFlag){print('no more edges remaining!')}
      break
    }
  }
  return(list(expandedV=expandedV, expandedE=expandedE, remainingE=remainingE))
}

IterateExpandGraph <- function(seedV, allE, maxK=Inf, printFlag=FALSE){
  out <- list()
  out$expandedV <- seedV
  out$expandedE <- NULL
  out$remainingE <- allE
  
  k <- 0
  while( (nrow(out$remainingE) > 0) && k<maxK ){
    
    out <- ExpandGraph(currentV=out$expandedV, currentE=out$expandedE, remainingE=out$remainingE)
    
    k <- k+1
    if(printFlag){print(sprintf('finished %d iterations', k))}
  }
  return(list(k=k, out=out))
}

# GetGeneReg <- function(edges=FALSE, toEntrez=TRUE){
#   W <- read.table(DataFile('tf_gene_network/wpg.txt'), header=TRUE, sep=' ')
#   
#   if(substr(colnames(W)[1],1,1) == 'X'){
#     colnames(W) <- sapply(colnames(W), function(x) substr(x, 2, nchar(x)))
#   }
#   
#   if(toEntrez){
#     uniprot <- rownames(W)
#     rownames(W) <- MapUniprot2Entrez(uniprot)
#   }
#   
#   if(edges){
#     out <- BinaryMatrix2Edges(W)
#   }else{
#     out <- W
#   }
#   return(out)
# }

GetChemStructure <- function(){
  #D <- read.table('../../data/chem_structure//lincs_feature.csv', sep=',', header=TRUE)
  #rownames(D) <- D[,1]
  #D <- D[,2:ncol(D)]
  n <- load(DataFile('chem_structure/chemStructure.RData'))
  stopifnot(n == 'chemStructure')
  return(chemStructure)
}

DataFile <- function(file){
  return(paste0(GetConfig('DATAPATH'),'/',file))
}

NumSigs <- function(T, dim='all'){
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

# GetTensorAnnot <- function(){
#    load(DataDir('/expr/tensor/tensor_annot.RData'))
#   return(annot)
# }

# GetDataTensor <- function(normalize=FALSE, dataset=NULL){
#   if(!is.null(dataset)){
#     load(DataDir(paste0('expr/tensor/', dataset, '.RData')))
#   }else{
#     load(ResultsDir('T_meas.RData'))
#     annot = GetTensorAnnot()
#     dimnames(T_meas) = list(drug=annot$pertIds, gene=annot$geneIds, cell=annot$cellIds)
#   }
# 
#   if(normalize){
#     T_meas = NormSigs(T_meas)
#   }
#   return(T_meas)
# }
# 
# GetAllTensors <- function(normalize=TRUE, meas=TRUE, cv=TRUE, comp=TRUE, merge=TRUE, pred=TRUE){
#   
#   if(all(c(!meas, !cv, !comp, !merge, !pred))){
#     warning('not returning any tensors!')
#   }
#   
#   tensors = list()
#   
#   X = GetDataTensor(normalize=FALSE)
#   dimnames = dimnames(X)
#   
#   if(meas){
#     tensors$meas = X
#   }
#   
#   for(typename in c('cv', 'comp', 'merge', 'pred')){
#     type = eval(parse(text=typename))
#     if(type==TRUE){
#       tensors[[typename]] = list()
#       for(method in c('mean', 'mean2', 'knn', 'tensor')){
#         load(sprintf('%s/T_%s_%s.RData', ResultsDir(), method, typename))
#         X = eval(parse(text=sprintf('T_%s_%s', method, typename)))
#         if(normalize){
#           X = NormSigs(X)
#         }
#         dimnames(X) = dimnames
#         tensors[[typename]][[method]] = X
#       }
#     }
#   }
#   
#   return(tensors)
# }

UnfoldTensor <- function(X, dim=1){
  return(wrap(X, map=list(dim,NA)))
}

# standardize tensor values per gene
TensorZ <- function(X, m=NULL, v=NULL){
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


TensorDEG <- function(X, normGene=FALSE, method='tensor', percDEG=2, symmetric=FALSE){
  if(all(c('drug', 'gene', 'cell') %in% names(dimnames(X)))){
    X = aperm(X, c('drug','gene', 'cell'))
  }

  if(normGene){
    # normalize each gene's values across tensor
    Z = TensorZ(X)$Z
  }else{
    # or not
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

CallDEG <- function(x, p, symmetric=FALSE){
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
#   print(sprintf('num DEGs called: %d out of %d', 
#                 length(which(abs(d)==1)), n_elts))
  return(d)
}


# FindSurprisingDEGs <- function(T_true, T_est, pertIds, geneIds, cellIds, n=100, maxPred = NULL, sortBy = 'geneId'){
#   a = abs(T_true) - abs(T_est)
#   b = sort(a, decreasing=TRUE)[1:n]
#   idx = which(a %in% b)
#   idx3 = arrayInd(idx,dim(a))
#   
#   out = list()
#   out$diff = a[idx]
#   out$true = T_true[idx]
#   out$pred = T_est[idx]
#   out$pertId = pertIds[idx3[,1]]
#   out$geneId = geneIds[idx3[,2]]
#   out$cellId = cellIds[idx3[,3]]
#   out = as.data.frame(out)
#   out = out[order(out[,sortBy]),]
#   
#   if(!is.null(maxPred)){
#     out = out[abs(out$pred) < maxPred,]
#   }
#   return(out)
# }


# DimStats <- function(A, dim, dimnames=NULL, plot=FALSE, orderBy=NULL){
#   B = wrap(A, map = list(dim, NA))
#   #B = abs(wrap(A, map = list(dim, NA)))
#   if(!is.null(dimnames)){
#     rownames(B) = dimnames
#   }else{
#     rownames(B) = dimnames(A)[[dim]]
#   }
# 
#   out = list()
#   out$mean = apply(B, 1, function(x) mean(x, na.rm=TRUE))
#   out$median = apply(B, 1, function(x) median(x, na.rm=TRUE))
#   out$min = apply(B, 1, function(x) min(x, na.rm=TRUE))
#   out$var = apply(B, 1, function(x) var(x, na.rm=TRUE))
#   out$se = apply(B, 1, function(x) std.error(x, na.rm=TRUE))
#   out$max1 = apply(B, 1, function(x) max(x, na.rm=TRUE))
#   out$max2 = apply(B, 1, function(x) sort(x, decreasing=TRUE)[2])
#   out$max3 = apply(B, 1, function(x) sort(x, decreasing=TRUE)[3])
#   out$max4 = apply(B, 1, function(x) sort(x, decreasing=TRUE)[4])
#   out$max5 = apply(B, 1, function(x) sort(x, decreasing=TRUE)[5])
#   
#   if(plot){
#     C = stack(as.data.frame(t(B)))
#     C$med = unlist(lapply(1:nrow(C), function(i) dim_med[C$ind[i]]))
#     C$ind2 = reorder(C$ind, C$med)
#     ggplot(C, aes(x = ind2, y = values)) +
#       geom_boxplot(fill='turquoise')
#   }
#   
#   out = as.data.frame(out)
#   
#   if(!is.null(orderBy)){
#     out = out[order(out[,orderBy]),]
#   }
# 
#   return(out)
# } 

ComputeAUC <- function(est, labels, computeROC=FALSE, abs=TRUE){
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

ComputeCellSpecificity <- function(X, normalize=FALSE){
  nDrug = dim(X)[1]
  nCell = dim(X)[3]
  # compute pairwise cosine distances between drug signatures in measured tensor
  A = array(data=NA, dim=c(nCell,nCell))
  cs = rep(NA, times=nDrug)
  r = list()
  for(d in 1:nDrug){
    print(d)
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

# GetCellSpecificity <- function(){
#   load(ResultsDir('cell_specificity_per_drug.RData'))
#   return(cellSpecificity)
# }

GetCancerDrivers <- function(){
  n = load('../../data/cancer_drivers/cancer_drivers_lm_intersect_v2.RData')
  stopifnot(n == 'drivers')
  return(drivers)
}

ComputeGeneGeneCor <- function(tensor, nGene=dim(tensor)[2], cellSpecific=FALSE){
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
      print(cell)
      out[[cell]] = cor(tensor[,1:nGene,cell], use='pairwise')
    }
  }
  return(out)
}

# GetGeneGeneCor <- function(){
#   load(ResultsDir('G.RData'))
#   return(G)
# }


RankSigs <- function(tensor){
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

# GetAccuracyPerMode <- function(){
#   load(ResultsDir('accuracy_per_mode.RData'))
#   return(C)
# }

# GetTensorDims <- function(dim=NA, fromFile=FALSE){
#   
#   if(fromFile){
#     X = GetDataTensor()
#     list[nDrug, nGene, nCell] = dim(X)
#   }else{
#     nDrug = 611
#     nGene = 978
#     nCell = 11
#   }
#   
#   if(is.na(dim)){
#     out = list(nDrug=nDrug, nGene=nGene, nCell=nCell)
#   }else{
#     if(dim == 'drug'){
#       out = nDrug
#     }else if(dim == 'gene'){
#       out = nGene
#     }else if(dim == 'cell'){
#       out = nCell
#     }else{
#       error('unexpected value for dim')
#     }
#   }
#   return(out)
# }

# GetExpCount <- function(){
#   load(DataDir('/expr/tensor/exp_count.RData'))
#   pertIds = GetTensorAnnot()$pertIds
#   exp_count = exp_count[pertIds,]
#   return(exp_count)
# }

CheckColumnStructure <- function(tensor){
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

RemoveMissingDataFromXY <- function(X, Y){
  idx1 = na.action(na.omit(X))
  idx2 = which(is.na(Y))
  idxNA = union(idx1, idx2)
  if(length(idxNA) > 0){
    X = as.matrix(X[-idxNA,])
    Y = Y[-idxNA]
  }
  return(list(X=X, Y=Y))
}

LoadTensorMat <- function(file){
  out = readMat(file)
  tensor = out$T
  geneIds = unlist(out$geneIds)
  pertIds = unlist(out$pertIds)
  cellIds = unlist(out$cellIds)
  dimnames(tensor) = list(pertIds, geneIds, cellIds)
  return(list(tensor=tensor, pertIds=pertIds, geneIds=geneIds, cellIds=cellIds))
}

GetGeneIdsTensor <- function(){
  load(DataDir('expr/tensor/d6_24hr_subsets/T50_1.RData'))
  return(dimnames(T_meas)[[2]])
}
