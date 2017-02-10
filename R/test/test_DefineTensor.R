TestSelectDataForTensor = function(){
  info = LoadCDInfo(debug=FALSE) # 
  pThresh = 0.0005
  print = FALSE
  removeDuplicates=FALSE
  
  # test individual filters and make sure the filters are satisfied
  out1 = SelectDataForTensor(info, pThresh=pThresh, print=print, removeDuplicates = removeDuplicates)
  stopifnot(max(out1$pvalue) <= pThresh)
  stopifnot(length(which(out1$pert_dose < 7)) > 0)
  stopifnot(length(which(out1$pert_time > 25)) > 0)
  stopifnot(length(unique(out1$cell_id)) > 10)
  stopifnot(length(unique(out1$pert_id)) > 1000)
  
  out2 = SelectDataForTensor(info, pThresh= pThresh, specificDose = TRUE, print=print, removeDuplicates = removeDuplicates)
  out2B = SelectDataForTensor(out1, pThresh= pThresh, specificDose = TRUE, print=print, removeDuplicates = removeDuplicates)
  stopifnot(identical(out2, out2B))
  stopifnot(nrow(out2) < nrow(out1))
  stopifnot(min(out2$pert_dose) >= 8 && max(out2$pert_dose) <= 12)


  out3 = SelectDataForTensor(info, pThresh=pThresh, specificDose=TRUE, time=6, print=print, removeDuplicates = removeDuplicates)
  out3B = SelectDataForTensor(out2, pThresh=pThresh, specificDose=TRUE, time=6, print=print, removeDuplicates = removeDuplicates)
  stopifnot(identical(out3, out3B))
  stopifnot(nrow(out3) < nrow(out2))
  stopifnot(sort(unique(out3$pert_time)) %in%  c(6, 6.4))

  out4 = SelectDataForTensor(info, pThresh=pThresh, specificDose=TRUE, time=6, nCells=5, print=print, removeDuplicates = removeDuplicates)
  out4B = SelectDataForTensor(out3, pThresh=pThresh, specificDose=TRUE, time=6, nCells=5, print=print, removeDuplicates = removeDuplicates)
  stopifnot(identical(out4, out4B))
  stopifnot(nrow(out4) < nrow(out3))
  stopifnot(length(unique(out4$cell_id)) == 5)
  
  out5 = SelectDataForTensor(info, pThresh=pThresh, specificDose=TRUE, time=6, nCells=5, nDrugs=20, print=print, removeDuplicates = removeDuplicates)
  out5B = SelectDataForTensor(out4, pThresh=pThresh, specificDose=TRUE, time=6, nCells=5, nDrugs = 20, print=print, removeDuplicates = removeDuplicates)
  stopifnot(identical(out5, out5B))
  stopifnot(nrow(out5) < nrow(out4))
  stopifnot(length(unique(out5$pert_id)) == 20)
  
  out6 = SelectDataForTensor(info, pThresh=pThresh, specificDose=TRUE, time=6, nCells=5, nDrugs=20, 
                             cellIds=c('MCF7', 'PC3', 'foo'), pertIds=c('BRD-A06352508', 'BRD-A17065207', 'bar'), print=print)
  stopifnot(length(unique(out6$pert_id)) == 2)
  stopifnot(length(unique(out6$cell_id)) == 2)
  
  out7 = SelectDataForTensor(info, pThresh=pThresh, specificDose=TRUE, time=6, cellIds=c('MCF7', 'PC3', 'foo'),
                             pertIds=c('BRD-A06352508', 'BRD-A17065207', 'bar'), print=print, removeDuplicates = removeDuplicates)
  stopifnot(identical(out6, out7))

  out8 = SelectDataForTensor(info, pThresh=pThresh, time='both', cellIds=c('MCF7', 'PC3', 'foo'),
                             pertIds=c('BRD-A06352508', 'BRD-A17065207', 'bar'), print=print, removeDuplicates = removeDuplicates)
  M = as.matrix(cast(out8[,c('pert_id', 'cell_id', 'pert_time')], pert_id ~ cell_id, function(x) length(unique(x)), value='pert_time'))
  stopifnot(all(sort(unique(as.numeric(M))) == c(0,2)))
  
  out9 = SelectDataForTensor(info, pThresh=pThresh, time='both', nPerDrug=5, nPerCell=10, print=FALSE, removeDuplicates = removeDuplicates)
  M = as.matrix(cast(out9[,c('pert_id', 'cell_id', 'pvalue')], pert_id ~ cell_id, function(x) length(x)>0, value='pvalue'))
  stopifnot(all(apply(M, 1, sum) >= 5))
  stopifnot(all(apply(M, 2, sum) >= 10))
  
} 

TestCountCellsPerDrug = function(){
  info = LoadCDInfo(debug=TRUE)
  sigs = LoadCDSigs(debug=TRUE)
  tensor = ConstructTensor(sigs=sigs, info=info, pThresh=1, specificDose=FALSE, time='all',
                        removeDuplicates=FALSE, debug=TRUE, print=FALSE)$tensor
  count1 = NumSigs(tensor,'drug')
  count2 = CountCellsPerDrug(subset(info, pvalue <= 1))[names(count1)]
  stopifnot(all(count1 == count2))
}

TestSummarizeMatrix = function(){
  # construct matrix with known output and test
  M = matrix(data=c(TRUE, FALSE, TRUE, FALSE, TRUE, TRUE), nrow=3, ncol=2)
  numSigs = 4
  numPossible = 6
  density = 2/3 * 100
  nDrugs = 3
  nCells = 2
  nPerDrug = 1
  nPerCell = 2

  out1 = SummarizeMatrix(M)
  out2 = sprintf('%d out of %d signatures (density %0.1f%%), covering %d drugs and %d cell lines. \nMin of %d sigs per drug and %d sigs per cell.',
                 numSigs, numPossible, density, nDrugs, nCells, nPerDrug, nPerCell)
  stopifnot(out1 == out2)
}



TestCompareTensors = function(){
  load(DataDir('tensors/T_test.RData'))
  T1 = T_test
  d = dim(T1)
  T2 = T1[sample(d[1]), sample(d[2]), sample(d[3])]
  stopifnot(CompareTensors(T1, T2))
  
  T1[1,1,1] = -T1[1,1,1]
  stopifnot(!CompareTensors(T1,T2))
  
  T1[1,1,1] = NA
  stopifnot(!CompareTensors(T1,T2))
}

TestConstructTensor = function(){
  # run ConstructTensor on small sample
  info = LoadCDInfo(debug=TRUE)
  sigs = LoadCDSigs(debug=TRUE)
  nDrugs = 5
  nCells = 2
  nGenes = 978
  pThresh = 0.1
  print=FALSE
  out = ConstructTensor(sigs=sigs, info=info, nDrugs=nDrugs, nCells=nCells, pThresh=pThresh, print=print, debug=TRUE)
  
  # test that all outputs have appropriate dimension
  stopifnot(all(dim(out$tensor) == c(nDrugs, nGenes, nCells)))
  stopifnot(all(dim(out$maxP) == c(nDrugs, nCells)))
  stopifnot(all(dim(out$nAvg) == c(nDrugs, nCells)))
  stopifnot(all(dim(out$nReps) == c(nDrugs, nCells)))
  stopifnot(length(out$nSigs) == 1)
  
  # test that all signatures are either norm 1 or all NA
  apply(out$tensor, c(1,3), function(x) stopifnot(all(is.na(x)) || abs(Norm2(x) - 1) < 1e-14))

  # test that all pvalues are <= pThresh and >= 0
  stopifnot(all(out$maxP <= pThresh, na.rm=TRUE))
  
  # test that nSigs equals the sum of nAvg
  stopifnot(out$nSigs == sum(out$nAvg))
  
  # test that nReps is always >= nAvg. And the sum is roughly greater by a factor of 3.
  stopifnot(all(out$nReps >= out$nAvg, na.rm=TRUE))
  stopifnot(mean(out$nReps / out$nAvg, na.rm=TRUE) > 2)
  
  # if I choose specificDose=TRUE, nAvg and nRep should be <= corresponding matrices from specificDose = FALSE
  out1 = ConstructTensor(sigs=sigs, info=info, nDrugs=nDrugs, nCells=nCells, pThresh=pThresh, specificDose=TRUE, print=print, debug=TRUE)
  out2 = ConstructTensor(sigs=sigs, info=info, nDrugs=nDrugs, nCells=nCells, pThresh=pThresh, specificDose=FALSE, print=print, debug=TRUE)
  stopifnot(all(out1$nAvg <= out2$nAvg, na.rm=TRUE))
  stopifnot(all(out1$nRep <= out2$nRep, na.rm=TRUE))
  
  # same for time = 6/24 vs. time='all'
  out1 = ConstructTensor(sigs=sigs, info=info, nDrugs=nDrugs, nCells=nCells, pThresh=pThresh, time=24, print=print, debug=TRUE)
  out2 = ConstructTensor(sigs=sigs, info=info, nDrugs=nDrugs, nCells=nCells, pThresh=pThresh, time = 'all', print=print, debug=TRUE)
  stopifnot(all(out1$nAvg <= out2$nAvg, na.rm=TRUE))
  stopifnot(all(out1$nRep <= out2$nRep, na.rm=TRUE))
  
  # test nPerDrug and nPerCell
  info = SelectDataForTensor(info=info, nDrugs=50, nCells=10, pThresh=pThresh, nPerDrug=2, nPerCell=3, print=print)
  out = ConstructTensor(sigs=sigs, info=info, nDrugs=50, nCells=10, pThresh=pThresh, nPerDrug=2, nPerCell=3, print=print, debug=TRUE)
  M = !is.na(out$tensor[,1,])
  stopifnot(min(apply(M, 1, sum)) >= 2)
  stopifnot(min(apply(M, 2, sum)) >= 3)
}

TestRestrictToCommonSigs = function(){
  # load three tensors where I know what the common support is
  load(DataDir('tensors/T_test.RData'))
  T1 = T_test[,,1:5]
  T2 = T_test[1:3,2:4,]
  T3 = T_test[3:1,1:4,]
  inList = list(T1=T1,T2=T2,T3=T3)
  outList = RestrictToCommonSigs(inList, print=FALSE)
  
  stopifnot(length(unique(outList))==1)
  stopifnot(all(dim(outList[[1]]) == c(3,3,5)))
  
  # now test a case where the differences lie in both missing signatures as well as different dimensions
  T4 = T_test + 2
  T4[1,,1] = NA # this one will remain in the outList
  T4[1,,6] = NA # we shouldn't see the effect of this one
  inList = list(T1=T1,T2=T2,T3=T3,T4=T4)
  outList2 = RestrictToCommonSigs(inList, print=FALSE)
  stopifnot(length(unique(outList2[1:3]))==1)
  stopifnot(length(unique(outList2))==2)
  stopifnot(NumSigs(outList2[[4]]) == NumSigs(outList[[1]]) - 1)
}

TestSubsetTensor = function(){
  load(DataDir('tensors/T_test.RData'))
  Tsm = SubsetTensor(T_test, nDrugs=4, nCells=7)
  stopifnot(CompareTensors(Tsm, T_test[c('BRD-A84481105','BRD-A52660433','BRD-A47513740','BRD-K55127134'),,
                                       c('A375','A549','HA1E','HT29','MCF7','PC3','VCAP')]))
}

TestSelectTopNDrugs = function(){
  info = LoadCDInfo(debug=TRUE)
  n = 100
  info2 = SelectTopNDrugs(info=info, nDrugs=n)
  stopifnot(length(unique(info2$pert_id)) == n)
  stopifnot(all(info[setdiff(rownames(info), rownames(info2)),'pert_id'] %ni% info2$pert_id))
}

TestSelectTopNCells = function(){
  info = LoadCDInfo(debug=TRUE)
  n = 1
  info2 = SelectTopNCells(info=info, nCells=n)
  stopifnot(length(unique(info2$cell_id)) == n)
  stopifnot(all(info[setdiff(rownames(info), rownames(info2)),'cell_id'] %ni% info2$cell_id))
}


TestLoadCDSigs = function(){
  # print(sprintf('Running first time..'))
  # sigs = LoadCDSigs(debug=FALSE) # takes about a minute to run
  # stopifnot(dim(sigs) == c(NDrugSigs(), 978)) 
  # sigs = LoadCDSigs(debug=FALSE) # This time it shouldn't reload, so should be almost instantaneous
  sigs = LoadCDSigs(debug=TRUE)
  stopifnot(dim(sigs) == c(1000,978))
}

TestLoadCDInfo = function(){
  info = LoadCDInfo(debug=TRUE)
  stopifnot(dim(info) == c(1000, 13))
  
  info = LoadCDInfo(debug=FALSE)
  stopifnot(dim(info) == c(NDrugSigs(), 13))
}

TestNDrugSigs = function(){
  stopifnot(is.numeric(NDrugSigs()))
}

TestEnsureNPerCell = function(){} #dummy
TestEnsureNPerDrug = function(){} #dummy
TestSummarizeInfo = function(){} #dummy
TestSummarizeTensor = function(){} #dummy
TestWriteTensor2Mat = function(){} #dummy
