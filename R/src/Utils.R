
#### directory and file handling ###################################################################

MakeDir = function(dir){
  if(!dir.exists(dir)){
    dir.create(dir, recursive=TRUE)
  }else{
    warning('directory already exists!')
  }
  stopifnot(dir.exists(dir))
}


PlotDir = function(file, subdir=NULL){
  stopifnot(identical(normalizePath(getwd()), normalizePath(BaseDir())))
  plotDir = paste('plot', DateStr(), subdir, sep='/')
  if(!file.exists(plotDir)){
    dir.create(plotDir, recursive=T)
  }
  return(paste(plotDir, file, sep='/'))
}

BaseDir = function(){
  return(GetConfig("LINCSPATH"))
}

GetConfig = function(key){
  a = system("cat config.txt;", intern=T)
  b = unlist(strsplit(a,"="))
  return(b[which(b==key)+1])
}

DateStr = function(){
  return(gsub('-', '_', Sys.Date()))
}

CheckDir = function(){
  dirs = unlist(strsplit(getwd(), "/"))
  stopifnot(identical(dirs[length(dirs)], 'dgc_predict'))
}

DataDir = function(subdir=''){
  dir = GetConfig('DATAPATH')
  return(paste0(dir, subdir))
}

ResultsDir = function(subdir=''){
  dir = GetConfig('RESULTSPATH')
  return(paste0(dir, subdir))
}

#### workspace helper functions ####################################################################

LsVars = function(envInt=0){
  setdiff(ls(env=sys.frame(envInt)), lsf.str(env=sys.frame(envInt)))
}

LsFcns = function(envInt=0){
  lsf.str(env=sys.frame(envInt))
}

GetVarName = function(var){
  return(deparse(substitute(var)))
}

#### set operations and helper functions ###########################################################

"%ni%" = Negate("%in%")

IsPartition = function(partitionList, fullSet){
  isPartition = TRUE
  
  for(i in 1:length(partitionList)){
    for(j in IncreasingSequence(i+1, length(partitionList))){
      if(IsIntersection(partitionList[[i]], partitionList[[j]])){
        isPartition = FALSE
      }
    }
  }
  
  all = partitionList[[1]]
  for(i in 2:length(partitionList)){
    all = union(all, partitionList[[i]])
  }
  
  if(!identical(as.numeric(sort(all)), as.numeric(sort(fullSet)))){
    isPartition = FALSE
  }
  return(isPartition)
}

IsIntersection = function(setA, setB){
  isIntersection = FALSE
  if(length(intersect(setA, setB))>0){
    isIntersection = TRUE
  }
  return(isIntersection)
}

# find all occurrences of each element of A in B
MatchAll = function(A, B, collapse='|'){
  out = c()
  for(a in A){
    matches = paste(which(B == a), collapse=collapse)
    out = c(out, matches)
  }
  out[out==''] = NA
  return(out)
}

# get elements that appear at least k times
GetMultiples = function(x, k=2){
  y = unique(x)
  out = c()
  for(a in y){
    count = length(which(x == a))
    if(count >= k){
      out = c(out, a)
    }
  }
  return(out)
}

# input is a list of items that may have some repeats
# output is a list of unique items, and the corresponding counts


#### math/stats stuff ##############################################################################

FisherExact = function(A, B, n_universe=length(union(A,B)), n_hypotheses=1, printFlag=T, alternative='greater'){
  aa = length(intersect(A, B))
  bb = length(setdiff(A, B))
  cc = length(setdiff(B, A))
  dd = n_universe - (aa + bb + cc)

  x = matrix(c(aa,bb,cc,dd), nrow=2, ncol=2, byrow=T)
  p = fisher.test(x, alternative=alternative)$p.value
  
  if(printFlag){
    print(sprintf('  %d in A, %d in B, %d in overlap, p=%f', length(A), length(B), aa, p))
  }
  return(list(p=p, adjp=min(p*n_hypotheses, 1), overlap=aa))
}

CosineDistance = function(x, y, normalize=TRUE, na.rm=TRUE){
  if(na.rm && normalize){
    idx1 = which(!is.na(x))
    idx2 = which(!is.na(y))
    idx = intersect(idx1,idx2)
    x = x[idx]
    y = y[idx]
    d = sum(x*y) / (sqrt(sum(x*x))*sqrt(sum(y*y)))
  }else{
    d = sum(x*y, na.rm=na.rm)
    if(normalize){
      d = d/(sqrt(sum(x*x))*sqrt(sum(y*y)))
    }
  }
  return(1-d)
}

# return a matrix with i,j entry corresponding to the cosine distance between rows i and j
CosineDistanceMatrix = function(M, normalize=TRUE, asVector=FALSE, fullMatrix=FALSE){
  D = M %*% t(M)
  
  if(normalize){
    norm = apply(M, 1, Norm2)
    N = norm %*% t(norm)
    D = D / N
  }
  
  if(!fullMatrix){
    for(i in 1:nrow(D)){
      for(j in 1:i){
        D[i,j] = NA
      }
    }
  }

  if(asVector){
    D = na.omit(as.vector(D))
  }
  
  if(!is.null(rownames(M))){
    rownames(D) = rownames(M)
    colnames(D) = rownames(M)
  }
  return(1-D)
}

RescaleVec = function(x, a=0, b=1, abs=TRUE){
  if(abs){
    x = abs(x)
  }
  r = range(x, na.rm=TRUE)
  y = a + (x - r[1])*(b - a) / (r[2] - r[1])
  return(y)
}

JaccardIndex = function(binary_vec1, binary_vec2){
  stopifnot(length(binary_vec1) == length(binary_vec2))
  
  idx1 = which(binary_vec1 == 1)
  idx2 = which(binary_vec2 == 1)
  
  return(JaccardIndex_fromIdx(idx1,idx2))
}

JaccardIndex_fromIdx = function(idx1, idx2){
  A_intersect_B = length(intersect(idx1, idx2))
  A_union_B = length(union(idx1, idx2))
  return(A_intersect_B / A_union_B)
}

ZScore = function(x){
  return( (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE))
}

#### data frame helper functions ###################################################################

FixNAStrings = function(df, printFlag=FALSE){
  for(i in 1:ncol(df)){
    idx = which(df[,i] == 'NA')
    df[idx,i] = NA
    if(printFlag){
      print(sprintf('fix na strings %d',i))
      print(idx)
    }
  }
  return(df)
}

ChangeColumnName = function(df, from, to){
  
  stopifnot(length(from) == length(to))
  if(length(from) > 1){
    for(i in 1:length(from)){
      df = ChangeColumnName(df, from[i], to[i])
    }
  }
  
  idx = which(names(df) == from)
  if(length(idx) != 1){
    warning(sprintf('found %d column name matches for %s', length(idx), from))
  }
  names(df)[idx] = to
  return(df)
}

RemoveDfColumns = function(df, columnNames){
  stopifnot(all(columnNames %in% names(df)))
  idx = which(names(df) %in% columnNames)
  return(df[,-idx, drop=F])
}

RemoveDfRows = function(df, rowNames){
  stopifnot(all(rowNames %in% rownames(df)))
  idx = which(rownames(df) %in% rowNames)
  return(df[-idx,])
}

Factor2Char = function(df){
  idx = sapply(df, is.factor)
  df[idx] = lapply(df[idx], as.character)
  return(df)
}

RemoveDuplicateRowsMulti = function(df, colNames){
  for(colName in colNames){
    df = RemoveDuplicateRows(df, colName)
  }
  return(df)
}

RemoveDuplicateRows = function(df, colName){
  x =  df[,colName]
  dups = unique(x[duplicated(x)])
  idxDups = which(x %in% dups)
  if(length(idxDups) > 0){
    out = df[-idxDups,]
  }else{
    out = df
  }
  return(out)
}

CollapseDuplicateRows = function(df, colName){
  x = df[,colName]
  idxDups = which(duplicated(x))
  if(length(idxDups) > 0){
    out = df[-idxDups,]
  }else{
    out = df
  }
  return(out)
}

MergeDfList = function(dfList, by, allRows=TRUE, allSuffixes=TRUE){
  if(allSuffixes){
    for(i in 1:length(dfList)){
      idx = which(names(dfList[[i]]) %in% by)
      names(dfList[[i]])[-idx] = paste(names(dfList[[i]])[-idx], names(dfList)[i], sep='.')
    }
  }
  mergedDf = as.data.frame(Reduce(function(x, y) merge(x, y, by=by, all=allRows), dfList))
  return(FixNAStrings(mergedDf))
}

CompareDfSubset = function(df1, df2, colNames=intersect(names(df1), names(df2)), matchName=NULL, 
                            printFlag=T, maxThresh=0, normThresh=0){
  
  stopifnot(all(colNames %in% names(df1)))
  stopifnot(all(colNames %in% names(df2)))
  
  namesA = setdiff(names(df1), names(df2))
  namesB = setdiff(names(df2), names(df1))
  
  df1 = df1[,colNames]
  df2 = df2[,colNames]
  
  if(!is.null(matchName)){
    idx = match(df1[,matchName], df2[,matchName])
    df2 = df2[idx,]
  }
  
  same = TRUE
  badNames = c()
  if(!identical(df1, df2)){
    for(name in colNames){
      if(!all(na.omit(df1[,name]) == na.omit(df2[,name]))){
        misMatches = which(na.omit(df1[,name]) != na.omit(df2[,name]))
        if(length(misMatches)>0){
          if(printFlag){print(name); cat('  mismatches: '); cat(misMatches); cat('\n')}
          if(!VectorCompare(df1[,name], df2[,name], maxThresh=maxThresh, normThresh=normThresh)){
            same = FALSE
            badNames = c(badNames, name)
          }
        }
      }
    }
  }
  
  goodNames = setdiff(colNames, badNames)
  
  return(list(same=same, badNames=badNames, goodNames=goodNames, namesA=namesA, namesB=namesB))
}

DataSummary = function(data, varname, groupnames){
  summary_func = function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))}
  data_sum = ddply(data, groupnames, .fun=summary_func, varname)
  data_sum = ChangeColumnName(data_sum, from='mean', to=varname)
  return(data_sum)
}

Write2XLS = function(dfnames, file, AdjWidth=TRUE, rownames=FALSE, colnames=TRUE){
  library(WriteXLS)
  WriteXLS(dfnames, row.names=rownames, col.names=colnames, AdjWidth=AdjWidth, FreezeCol=0, FreezeRow=1, BoldHeaderRow = TRUE,
           ExcelFileName=file)
}

ReadXLS = function(filename, sheetname=NULL) {
  sheets = readxl::excel_sheets(filename)
  x = lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  names(x) = sheets
  if(!is.null(sheetname)){
    if(sheetname %in% sheets){
      out = x[[sheetname]]
    }else{
      warning('sheetname not found; returning all sheets.')
    }
  }else{
    out = x
  }
  
  if(length(out) == 1){
    out = out[[1]]
  }
  return(out)
}

#### miscellaneous #################################################################################

Eval = function(string){
  return(eval(parse(text=string)))
}

AlphaNames = function(n){
  stopifnot(n <= 26^2)
  return(as.character(levels(interaction(letters, LETTERS)))[1:n])
}

Str2Vec = function(strings, split=',', unique=FALSE){
  if(unique){
    out = lapply(strings, function(x) unique(unlist(strsplit(x, split=split))))
  }else{
    out = lapply(strings, function(x) unlist(strsplit(x, split=split)))
  }
  if(length(strings)==1){
    out = unlist(out)
  }
  if(!is.null(names(strings))){
    names(out) = names(strings)
  }
  return(out)
}

Ind2Sub = function(num.columns, ind){
  c = ((ind-1) %% num.columns) + 1
  r = floor((ind-1) / num.columns) + 1
  return(cbind(r,c))
}

IncreasingSequence = function(from, to){
  if (to >= from){
    out = from:to
  }else{
    out = seq(1,1,length=0)
  }
  return(out)
}

GetVarianceQuantile = function(M, quantile){
  v = rowVars(M)
  idx = which(v > quantile(v, probs=(1-quantile)))
  names(idx) = NULL
  return(idx)
}

GetVarianceTopK = function(M, k){
  v = rowVars(M)
  return(which(rank(-v) <= k))
}

GetMedianTopK = function(M, k){
  m = rowMedians(M)
  return(which(rank(-m) <= k))
}

IsSymmetric = function(A, thresh=1e-12){
  symmetric = TRUE
  if(nrow(A) != ncol(A)){
    symmetric = FALSE
  }
  a = A - t(A)
  if(norm(a) > thresh){
    symmetric = FALSE
  }
  return(symmetric)
}

#### vector manipulation functions #################################

Norm2 = function(x, na.rm=TRUE){
  return(sqrt(sum(x^2, na.rm=na.rm)))
}

VectorCompare = function(x, y, normThresh=1e-12, maxThresh=NULL){
  isSame = TRUE
  if((length(x) != length(y)) || (is.na(x) != is.na(y)) || 
       class(x) != class(y) || any(is.finite(x) != is.finite(y)) ){
    isSame = FALSE
  }
  
  x = na.omit(x)
  y = na.omit(y)
  
  x = x[is.finite(x)]
  y = y[is.finite(y)]
  
  if(length(x) != length(y)){
    isSame = FALSE
  }
  
  if(isSame && is.numeric(x) && is.numeric(y) && length(x) > 0){
    
    if(is.null(normThresh) && is.null(maxThresh)){
      stop('both normThresh and maxThresh can\'t be null')
    }
    
    if(!is.null(normThresh)){
      if(Norm2(x - y) > normThresh){
        isSame = FALSE
      }
    }
    
    if(!is.null(maxThresh)){
      if(max(abs(x - y)) > maxThresh){
        isSame = FALSE
      }
    }
    
  }else if(isSame && !all(x == y)){
    isSame = FALSE
  }
  
  return(isSame)
}

FindAllDuplicates = function(v){
  return(which(duplicated(v) | duplicated(v, fromLast=TRUE)))
}

#### matrix manipulation functions #################################
MatrixCast = function(M, type){
  out = t(apply(M, 1, eval(parse(text=paste0('as.', type)))))
  stopifnot(dim(out) == dim(M))
  rownames(out) = rownames(M)
  colnames(out) = colnames(M)
  return(out)
}

SymmetrifyMatrix = function(M, diag=NA){
  if(nrow(M) != ncol(M)){
    warning('matrix is not symmetric, cannot symmetrify')
    out = M
  }else{
    M[is.na(M)] = 0
    M_t = t(M)
    diag(M_t) = 0
    out = M + M_t
    if(is.numeric(diag)){
      diag(out) = diag
    }
  }
  return(out)
}

GetUpperTriVec = function(M){
  return(as.vector(M[upper.tri(M, diag = FALSE)]))
}

MatrixCor = function(M1, M2, symmetric=TRUE){
  stopifnot(dim(M1) == dim(M2))
  if(symmetric){
    M1 = GetUpperTriVec(M1)
    M2 = GetUpperTriVec(M2)
  }
  return(cor(as.vector(M1), as.vector(M2), use='pairwise'))
}

CorMatrixList = function(list1, list2){
  stopifnot(length(list1) == length(list2))
  out = sapply(1:length(list1), function(i) MatrixCor(list1[[i]], list2[[i]]))
  names(out) = names(list1)
  return(out)
}

ClusterRows = function(M){
  hcr = hclust(dist(M))
  ddr = as.dendrogram(hcr)
  Rowv = rowMeans(M, na.rm=TRUE)
  ddr = reorder(ddr, Rowv)
  rowInd = order.dendrogram(ddr)
  return(rowInd)
}

ClusterCols = function(M){
  return(ClusterRows(t(M)))
}

ClusterRowsAndCols = function(M){
  return(list(rowInd=ClusterRows(M), colInd=ClusterCols(M)))
}


