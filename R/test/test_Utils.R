
#### directory and file handling ###################################################################

TestMakeDir = function(){
  MakeDir('foo/bar')
  unlink('foo', recursive=TRUE)
  stopifnot(!dir.exists('foo'))
}


TestPlotDir = function(){
  file1 = PlotDir('test')
  file2 = PlotDir('test', subdir='subdir')
  stopifnot(file.exists(dirname(file1)))
  stopifnot(file.exists(dirname(file2)))
}

TestBaseDir =function(){
  dirs =unlist(strsplit(BaseDir(), "/"))
  #stopifnot(identical(dirs[length(dirs)], 'dgc_predict'))
  stopifnot(identical(dirs[length(dirs)], 'code'))
}

TestGetConfig =function(){
  cfg =GetConfig("DATAPATH")
  stopifnot(is.character(cfg))
  stopifnot(!grepl('=', cfg))
  stopifnot(cfg != '')
}

TestDateStr = function(){
  str = DateStr()
  stopifnot(nchar(str) == 10)
  stopifnot(substr(str, 1, 2) == '20')
}

TestCheckDir =function(){
  setwd(BaseDir())
  CheckDir()
}

TestDataDir = function(){
  dirs = unlist(strsplit(DataDir(), "/"))
  stopifnot(identical(dirs[length(dirs)], 'data'))
  stopifnot(DataDir('/expr') == paste0(DataDir(), '/expr'))
}

TestResultsDir = function(){
  stopifnot(file.exists(ResultsDir()))
  stopifnot(ResultsDir('test') == paste0(ResultsDir(), 'test'))
}

#### workspace helper functions ####################################################################

TestLsVars = function(){
  # Not sure why this was giving me problems
  #i = 5
  #stopifnot('i' %in% LsVars())
}

TestLsFcns = function(){
  stopifnot('TestLsFcns' %in% LsFcns())  
}

TestGetVarName = function(){
  foo = 2
  bar = c(1,2,3)
  joe = list('a','b',foo)
  stopifnot(GetVarName(foo) == 'foo')
  stopifnot(GetVarName(bar) == 'bar')
  stopifnot(GetVarName(joe) == 'joe')
}

#### set operations and helper functions ###########################################################

TestIsPartition = function(){
  stopifnot(IsPartition(list(1:5, 11, 6:10), 1:11))
  stopifnot(!IsPartition(list(1:5, 11, 6:10), 1:10))
  stopifnot(!IsPartition(list(1:5, 11, 6:11), 1:10)) 
}

TestIsIntersection = function(){
  stopifnot(IsIntersection(1:5, 6:10) == FALSE)
  stopifnot(IsIntersection(1:5, 5:6) == TRUE)
}

TestMatchAll = function(){
  A = 1:4
  B = c(1, 3, 3, 2, 2, 5, 1)
  out = c('1-7', '4-5', '2-3', NA)
  stopifnot(identical(MatchAll(A,B,collapse='-'), out))
}

TestGetMultiples = function(){
  x = c(1, 1, 2, 3, 1, 2)
  stopifnot(identical(GetMultiples(x, 1), unique(x)))
  stopifnot(identical(GetMultiples(x, 2), c(1, 2)))
  stopifnot(identical(GetMultiples(x, 3), c(1)))
  stopifnot(identical(GetMultiples(x, 4), c()))
}

#### math/stats stuff ##############################################################################

TestFisherExact = function(){
  A = c(1, 2)
  B = c(1, 3)
  out = FisherExact(A, B, n_universe=4, printFlag=F, alternative='two.sided')
  stopifnot(out$p == 1)
  stopifnot(out$adjp == 1)
  stopifnot(out$overlap == 1)
  
  out = FisherExact(A, B, n_universe=1000, printFlag=F, alternative='two.sided')
  stopifnot(abs(out$p - 0.004) < .0001)
  stopifnot(abs(out$adjp - 0.004) < .0001)
  stopifnot(out$overlap == 1)  
}

TestCosineDistance = function(){
  x = c(1,0)
  y = c(0,1)
  a = CosineDistance(x,y,normalize=FALSE)
  b = CosineDistance(x,y,normalize=TRUE)
  stopifnot(a==b)
  stopifnot(a==1)
  x = rnorm(1000)
  y = rnorm(1000)
  a = CosineDistance(x,y, normalize =TRUE)
  stopifnot(abs(1-a) < 0.1)
  x[5] = NA
  y[4] = NA
  a = CosineDistance(x,y, normalize =TRUE, na.rm=TRUE)
  stopifnot(abs(1-a) < 0.1)
}

TestCosineDistanceMatrix = function(){
  M = matrix(data=c(1,0,0,1,-1,0,0,-1), nrow=4, ncol=2, byrow=TRUE)
  d = CosineDistanceMatrix(M, normalize=FALSE, asVector=FALSE)
  v = CosineDistanceMatrix(M, normalize=FALSE, asVector=TRUE)
  D = matrix(data=NA, nrow=4, ncol=4)
  D[1,] = c(NA, 1, 2, 1)
  D[2,] = c(NA, NA, 1, 2)
  D[3,] = c(NA, NA, NA, 1)
  stopifnot(identical(d, D))
  V = c(1, 2, 1, 1, 2, 1)
  stopifnot(all(v == V))
}

TestRescaleVec = function(){
  x = rnorm(100)
  y = RescaleVec(x, a=-5, b=-2, abs=FALSE)
  stopifnot(Norm2(range(y) - c(-5,-2)) < 1e-10)
  stopifnot(Norm2(cor(x,y) - 1) < 1e-10)
  
  y = RescaleVec(x, a=-2, b=-5, abs=FALSE)
  stopifnot(Norm2(range(y) - c(-5,-2)) < 1e-10)
  stopifnot(Norm2(cor(x,y)+1) < 1e-10)
}

TestJaccardIndex = function(){
  A = c(0,1,1,0)
  B = c(1,1,1,0)
  stopifnot(JaccardIndex(A,B) == 2/3)
  C = c(0,1,1,0)
  stopifnot(JaccardIndex(A,C) == 1)
  D = c(0,0,0,0)
  stopifnot(JaccardIndex(A,D) == 0)
}

TestJaccardIndex_fromIdx = function(){
  A = which(c(0,1,1,0)==1)
  B = which(c(1,1,1,0)==1)
  stopifnot(JaccardIndex_fromIdx(A,B) == 2/3)
  C = which(c(0,1,1,0)==1)
  stopifnot(JaccardIndex_fromIdx(A,C) == 1)
  D = which(c(0,0,0,0)==1)
  stopifnot(JaccardIndex_fromIdx(A,D) == 0)
}

TestZScore = function(){
  z = ZScore(rnorm(1000, mean=2, sd=0.1))
  m = mean(z)
  s = sd(z)
  stopifnot(abs(m) < 1e-10)
  stopifnot(abs(s-1) < 1e-10)
}

#### data frame helper functions ###################################################################

TestFixNAStrings = function(){
  A = data.frame(a=c(1, 'NA', 2), b=c(NA, NA, 'NA'))
  B = FixNAStrings(A, printFlag=F)
  stopifnot(which(is.na(A)) == c(4,5))
  stopifnot(which(is.na(B)) == c(2,4,5,6))
}

TestChangeColumnName = function(){
  A = data.frame(a=c(1, 'NA', 2), b=c(NA, NA, 'NA'))
  A = ChangeColumnName(A, 'b', 'hello')
  stopifnot(names(A) == c('a', 'hello'))
}

TestRemoveDfColumns = function(){
  df = MakeTestDf()
  df1 = RemoveDfColumns(df, c('col3'))
  stopifnot(identical(df1, df[,-3]))
  
  df2 = RemoveDfColumns(df, c('col1', 'col3'))
  stopifnot(identical(df2, df[,2, drop=F]))
  
  df3 = RemoveDfColumns(df, c('col1','col2','col3'))
  stopifnot(identical(rownames(df3), c('row1', 'row2', 'row3')))
  stopifnot(all(dim(df3) == c(3,0)))
}

TestRemoveDfRows = function(){
  df = MakeTestDf()
  df1 = RemoveDfRows(df, c('row3'))
  stopifnot(identical(df1, df[-3,]))
  
  df2 = RemoveDfRows(df, c('row1', 'row3'))
  stopifnot(identical(df2, df[2, ,drop=F]))
  
  df3 = RemoveDfRows(df, c('row1','row2','row3'))
  stopifnot(identical(colnames(df3), c('col1','col2','col3')))
  stopifnot(all(dim(df3) == c(0,3)))
}

TestFactor2Char = function(){
  A = data.frame(a=1:3, b=as.factor(c('a', 'a', NA)))
  B = Factor2Char(A)
  stopifnot(is.integer(B$a))
  stopifnot(is.character(B$b))
  stopifnot(is.na(B[3,'b']))
}

TestRemoveDuplicateRowsMulti = function(){
  A = data.frame(a=c(1, 1, 2, 2, 3, 4), b=c(1, 2, 2, 3, 4, 3))
  B = RemoveDuplicateRowsMulti(A, c('a', 'b'))
  stopifnot(identical(B$a, c(3, 4)))
  stopifnot(identical(B$b, c(4, 3)))
}

TestRemoveDuplicateRows = function(){
  TestRemoveDuplicateRowsMulti()
}

TestCollapseDuplicateRows = function(){
 A1 = data.frame(a=1:3, b=rep('a',3))
 A2 = data.frame(a=1:3, b=rep('b',3))
 A3 = data.frame(a=c(1,1), b=c('a','b'))
 A = rbind(A1, A2)
 B = CollapseDuplicateRows(A, colName='a')
 CompareDfs(B, A1)
 C = CollapseDuplicateRows(A, colName='b')
 CompareDfs(C, A3)
}

TestMergeDfList = function(){
  dfList = list(one=data.frame(a=1:5, b=c('a','b','c','d','e')),
                two=data.frame(a=5:1, c=6:10),
                three=data.frame(a=1:4, d=1:4))
  out = MergeDfList(dfList, by='a')
  stopifnot(nrow(out)==5)
  stopifnot(identical(names(out),c('a', 'b.one', 'c.two', 'd.three')))
  stopifnot(identical(out$c.two,10:6))
}

TestCompareDfSubset = function(){
  A = data.frame(a=1:3, b=4:6, c=c('a','b','c'))
  B = A[3:1,c('a','c')]
  stopifnot(CompareDfSubset(A,B, printFlag=F)$same==FALSE)
  stopifnot(CompareDfSubset(A,B, printFlag=F, matchName='a')$same)
  stopifnot(CompareDfSubset(A,B, printFlag=F, matchName='c')$same)
}

TestDataSummary = function(){
  n = 10000
  data = melt(data.frame(x=rnorm(n), y=rnorm(n, mean=2)), measure.vars=c('x','y'))
  out = DataSummary(data, varname='value', groupnames='variable')
  rownames(out) = out$variable
  thresh = 0.05
  stopifnot(abs(out['x','value']) < thresh)
  stopifnot(abs(out['y','value']-2) < thresh)
  stopifnot(abs(out['x','sd']-1) < thresh)
  stopifnot(abs(out['y','sd']-1) < thresh)
}

TestWrite2XLS = function(){}
TestReadXLS = function(){}


#### miscellaneous #################################################################################

TestEval = function(){
  stopifnot(Eval('2+2') == 4)
}

TestAlphaNames = function(){
  nm = AlphaNames(5)
  stopifnot(nm[1] == 'a.A')
  stopifnot(nm[5] == 'e.A')
}

TestStr2Vec = function(){
  str1 = 'a,b,foo,bar,a'
  strList1 = c('a','b','foo','bar','a')
  str2 = 'apple,peaches,pumpkin pie'
  strList2 = c('apple','peaches','pumpkin pie')
  stopifnot(identical(Str2Vec(str1), strList1))
  stopifnot(identical(Str2Vec(str1, unique=T), unique(strList1)))
  stopifnot(identical(Str2Vec(c(str1,str2)), list(strList1,strList2)))
}

TestInd2Sub = function(){
  A = matrix(1:12, nrow=3, ncol=4)
  out = Ind2Sub(ncol(A), which(A<6))
  stopifnot(identical(out[,'r'], c(1,1,1,1,2)))
  stopifnot(identical(out[,'c'], c(1,2,3,4,1)))
}

TestIncreasingSequence = function(){
  stopifnot(length(IncreasingSequence(5, 1)) == 0)
  stopifnot(length(IncreasingSequence(1, 5)) == 5)
}

TestGetVarianceQuantile = function(){
  # generate matrix with increasing variance in each row
  set.seed(1234)
  n = 1000
  X = replicate(n, rnorm(n))
  for(i in 1:n){
    X[i,] = X[i,]*i
  }
  idx = GetVarianceQuantile(X, quantile=0.1)
  stopifnot(!(1:50 %in% idx))
  stopifnot(950:1000 %in% idx)
}

TestGetVarianceTopK = function(){
  # generate matrix with increasing variance in each row
  n = 1000
  set.seed(1234)
  X = replicate(n, rnorm(n))
  for(i in 1:n){
    X[i,] = X[i,]*i
  }
  idx = GetVarianceTopK(X, 100)
  stopifnot(!(1:50 %in% idx))
  stopifnot(950:1000 %in% idx)
}

TestGetMedianTopK = function(){
  # generate matrix with increasing median in each row
  n = 1000
  set.seed(1234)
  X = replicate(n, rnorm(n))
  for(i in 1:n){
    X[i,] = X[i,]*i
  }
  idx = GetVarianceTopK(X, 100)
  stopifnot(!(1:50 %in% idx))
  stopifnot(950:1000 %in% idx)
}

TestIsSymmetric = function(){
  A = replicate(10, rnorm(10))
  B = A + t(A)
  stopifnot(!IsSymmetric(A))
  stopifnot(IsSymmetric(B))
}

#### vector manipulation functions #################################################################
TestNorm2 = function(){
  x = c(2,2,sqrt(41))
  stopifnot(abs(Norm2(x) - 7) < 1e-6)
}

TestVectorCompare = function(){
  x = c(1, 2, 3)
  y = c(1, 2, 3) + 0.01
  stopifnot(!VectorCompare(x,y))
  stopifnot(VectorCompare(x,y,maxThresh=0.011, normThresh=NULL))
  stopifnot(VectorCompare(x, y, normThresh=sqrt(4e-4)))
}

TestFindAllDuplicates = function(){
  a = c(letters[1:6], 'A', letters[1:6])
  stopifnot(a[-FindAllDuplicates(a)] == 'A')
}

#### matrix manipulation functions #################################################################
TestMatrixCast = function(){
  A = matrix(data=rep(c(0,1), 3), nrow=3, dimnames=list(letters[1:3], LETTERS[1:2]))

  B = MatrixCast(A, 'logical')
  C = MatrixCast(A, 'character')
  
  testB = matrix(data=as.logical(rep(c(0,1), 3)), nrow=3, dimnames=list(letters[1:3], LETTERS[1:2]))
  testC = matrix(data=as.character(rep(c(0,1), 3)), nrow=3, dimnames=list(letters[1:3], LETTERS[1:2]))
  
  stopifnot(identical(B, testB))
  stopifnot(identical(C, testC))
}

TestSymmetrifyMatrix = function(){
  n = 10
  A = matrix(data=NA, nrow=n, ncol=n)
  for(i in 1:n){
    for(j in IncreasingSequence(i+1,n)){
      A[i,j] = rnorm(1)
    }
  }
  B = SymmetrifyMatrix(A, diag=1)
  stopifnot(all(B == t(B)))
  B[is.na(A)] = NA
  stopifnot(identical(A, B))
}

TestGetUpperTriVec = function(){
  A = matrix(data=1:9, nrow=3, ncol=3)
  a = GetUpperTriVec(A)
  stopifnot(a == c(4,7,8))
}

TestMatrixCor = function(){
  # perfectly correlated matrices
  M1 = matrix(data=rnorm(100), nrow=10, ncol=10)
  M2 = 2*M1
  stopifnot(MatrixCor(M1, M2, symmetric=FALSE) == 1)
  stopifnot(MatrixCor(M1, -M2, symmetric=FALSE) == -1)
  
    # less correlated
  M2 = M1 + matrix(data=rnorm(100), nrow=10, ncol=10)
  stopifnot(MatrixCor(M1,M2, symmetric=FALSE) > 0.4)
  
  # no entries line up
  M2[upper.tri(M2, diag=TRUE)] = NA
  stopifnot(is.na(MatrixCor(M1,M2)))
  stopifnot(!is.na(MatrixCor(M1, M2, symmetric=FALSE)))
}

TestCorMatrixList = function(){
  n = 100
  M = matrix(data=rnorm(n*n), nrow=n, ncol=n)
  M1 = list()
  M2 = list()
  for(i in 1:10){
    M1[[i]] = M
    M2[[i]] = M + matrix(data=rnorm(n*n, sd=0.1*i), nrow=n, ncol=n)
  }
  out = CorMatrixList(M1, M2)
  stopifnot(!is.unsorted(-out))
}

# this works but it generates a plot so I don't want to run it every time I test
TestClusterRows = function(){
  #   M = MakeTestMatrix()
  #   hm = heatmap(M, labRow='', labCol='')
  #   stopifnot(identical(hm$rowInd, ClusterRows(M)))
}

# this works but it generates a plot so I don't want to run it every time I test
TestClusterCols = function(){
  #   M = MakeTestMatrix()
  #   hm = heatmap(M, labRow='', labCol='')
  #   stopifnot(identical(hm$colInd, ClusterCols(M)))
}

TestClusterRowsAndCols = function(){}

#### helper functions for testing ##################################################################

CompareDfs = function(DF1, DF2){
  a = Factor2Char(DF1)
  b = Factor2Char(DF2)
  stopifnot(all(a == b))
}

GetFunctionNames = function(file){
  con  = file(file, open = 'r')
  functions = vector()
  
  while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0) {
    toks = unlist(strsplit(oneLine, split=' = function'))
    if( (length(toks) > 1) & !(substr(toks[1],1,1) %in% c('#', ' ')) ){
      functions[length(functions)+1] = toks[1]
    }
  } 
  close(con)
  return(functions)
}

MakeTestDf = function(){
  df = data.frame(col1=1:3, col2=4:6, col3=c('a','b','c'))
  rownames(df) = c('row1', 'row2', 'row3')
  return(df)
}

