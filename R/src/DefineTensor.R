
# filter CD signatures to be input into the tensor
# info should be the data frame stored in DataDir('expr/rdata/allDrugSigs.RData')

SelectDataForTensor = function(info, pThresh=1, specificDose=FALSE,
                                time='all', nCells=NA, nDrugs=NA, 
                                cellIds = NA, pertIds = NA, 
                                nPerDrug=1, nPerCell=1, print=TRUE,
                                removeDuplicates=TRUE, annot = GetLincsAnnot()){
  
  # filter by p-value
  info = subset(info, pvalue <= pThresh)
  if(print){print(sprintf('After filtering by pvalue: %s', SummarizeInfo(info)))}
  
  # filter by dose
  if(specificDose){
    info = subset(info, pert_dose > 8 & pert_dose < 12)
    if(print){print(sprintf('After filtering by dose: %s', SummarizeInfo(info)))}
  }
  
  # filter by time point
  if(time == 6 || time == 24){
    info = subset(info, round(pert_time) == time)
    if(print){print(sprintf('After filtering by time point: %s', SummarizeInfo(info)))}
  }else if (time == 'both'){
    info = subset(info, round(pert_time) %in% c(6,24))
    
    # new columns
    pert_cell = paste0(info$pert_id, info$cell_id)
    rounded_time = round(info$pert_time)
    
    # create a table to count unique timepoints
    counts = aggregate(rounded_time, by=list(pert_cell), function(x) length(unique(x)))
    names(counts) = c('pert_cell', 'time_count')
    
    # identify list of unique drug/cells to keep
    counts = subset(counts, time_count==2)
    
    # subset based on this list
    info = subset(info, pert_cell %in% counts$pert_cell)
    if(print){print(sprintf('After filtering by time point: %s', SummarizeInfo(info)))}
  }else if (time != 'all'){
    stop('unrecognized input')
  }
  
  if(any(!is.na(cellIds))){
    info = subset(info, cell_id %in% cellIds)
    if(print){print(sprintf('After filtering to chosen cell ids: %s', SummarizeInfo(info)))}
  }
  
  if(any(!is.na(pertIds))){
    info = subset(info, pert_id %in% pertIds)
    if(print){print(sprintf('After filtering to chosen pert ids: %s', SummarizeInfo(info)))}
  }
  
  #### Remove duplicate perts
  if(removeDuplicates){
    uniquePerts = unique(info$pert_id)
    drugNames = annot$name[match(uniquePerts, annot$pert_id)]
    countPerDrug = CountCellsPerDrug(info)[uniquePerts]
    pertDf = data.frame(pert_id=uniquePerts, name=drugNames, n_cells=countPerDrug)
    pertDfKeep = ddply(pertDf,  'name', summarise, pert_id=pert_id[which.max(n_cells)])
    info = subset(info, pert_id %in% pertDfKeep$pert_id)
  }
 
  #### FIRST TRY ONE ORDERING
  infoSave = info
  
  info = EnsureNPerCell(info, nPerCell)
  if(print){print(sprintf('After filtering by nPerCell: %s', SummarizeInfo(info)))}
  
  info = EnsureNPerDrug(info, nPerDrug)
  if(print){print(sprintf('After filtering by nPerDrug: %s', SummarizeInfo(info)))}
  
  info = SelectTopNCells(info, nCells)
  if(print){print(sprintf('After filtering by nCells: %s', SummarizeInfo(info)))}
  
  info = SelectTopNDrugs(info, nDrugs)
  if(print){print(sprintf('After filtering by nDrugs: %s', SummarizeInfo(info)))}

  # Check that SelectTopNCell/Drug functions didn't ruin NPerCell/Drug properties
  while(!identical(info, EnsureNPerCell(info, nPerCell)) || !identical(info, EnsureNPerDrug(info, nPerDrug))){
    info = EnsureNPerCell(info, nPerCell)
    if(print){print(sprintf('After filtering by nPerCell again: %s', SummarizeInfo(info)))}
    info = EnsureNPerDrug(info, nPerDrug)
    if(print){print(sprintf('After filtering by nPerDrug again: %s', SummarizeInfo(info)))}
  }
  
  info1 = info
  
  #### NOW TRY SECOND ORDERING
  
  info = SelectTopNCells(infoSave, nCells)
  if(print){print(sprintf('After filtering by nCells: %s', SummarizeInfo(info)))}
  
  info = SelectTopNDrugs(info, nDrugs)
  if(print){print(sprintf('After filtering by nDrugs: %s', SummarizeInfo(info)))}
  
  info = EnsureNPerCell(info, nPerCell)
  if(print){print(sprintf('After filtering by nPerCell: %s', SummarizeInfo(info)))}
  
  info = EnsureNPerDrug(info, nPerDrug)
  if(print){print(sprintf('After filtering by nPerDrug: %s', SummarizeInfo(info)))}
  
  # Check that SelectTopNCell/Drug functions didn't ruin NPerCell/Drug properties
  while(!identical(info, EnsureNPerCell(info, nPerCell)) || !identical(info, EnsureNPerDrug(info, nPerDrug))){
    info = EnsureNPerCell(info, nPerCell)
    if(print){print(sprintf('After filtering by nPerCell again: %s', SummarizeInfo(info)))}
    info = EnsureNPerDrug(info, nPerDrug)
    if(print){print(sprintf('After filtering by nPerDrug again: %s', SummarizeInfo(info)))}
  }
  
  info2 = info
  
  #### NOW COMPARE WHICH ONE GIVES ME MORE DATA
  if(!identical(info1, info2)){
    stop('Results from two data selection orderings is not the same, please check.')
  }
  
  return(info1)
}

SelectTopNCells = function(info, nCells){
  if(!is.na(nCells)){
    infoPerCell = split(info, f=info$cell_id)
    countPerCell = sapply(1:length(infoPerCell), function(i) length(unique(infoPerCell[[i]]$pert_id)))
    nCells = min(length(countPerCell), nCells)
    names(countPerCell) = names(infoPerCell)
    cells = names(sort(countPerCell, decreasing=TRUE))[1:nCells]
    info = subset(info, cell_id %in% cells)
  }
  return(info)
}

SelectTopNDrugs = function(info, nDrugs){
  if(!is.na(nDrugs)){
    infoPerDrug = split(info, f=info$pert_id)
    countPerDrug = sapply(1:length(infoPerDrug), function(i) length(unique(infoPerDrug[[i]]$cell_id)))
    nDrugs = min(length(countPerDrug), nDrugs)
    names(countPerDrug) = names(infoPerDrug)
    drugs = names(sort(countPerDrug, decreasing=TRUE))[1:nDrugs]
    info = subset(info, pert_id %in% drugs)
  }
  return(info)
}

EnsureNPerCell = function(info, nPerCell){
  if(!is.na(nPerCell)){
    infoPerCell = split(info, f=info$cell_id)
    countPerCell = sapply(1:length(infoPerCell), function(i) length(unique(infoPerCell[[i]]$pert_id)))
    names(countPerCell) = names(infoPerCell)
    cells = names(which(countPerCell >= nPerCell))
    info = subset(info, cell_id %in% cells)
  }
  return(info)
}

EnsureNPerDrug = function(info, nPerDrug){
  if(!is.na(nPerDrug)){
    countPerDrug = CountCellsPerDrug(info)
    drugs = names(which(countPerDrug >= nPerDrug))
    info = subset(info, pert_id %in% drugs)
  }
  return(info)
}

CountCellsPerDrug = function(info){
  infoPerDrug = split(info, f=info$pert_id)
  countPerDrug = sapply(1:length(infoPerDrug), function(i) length(unique(infoPerDrug[[i]]$cell_id)))
  names(countPerDrug) = names(infoPerDrug)
  return(countPerDrug)
}

SummarizeInfo = function(info){
  return(sprintf('%d signatures, %d unique drugs, %d unique cell types', nrow(info), 
                 length(unique(info$pert_id)), length(unique(info$cell_id))))
}

SummarizeMatrix = function(M){
  stopifnot(is.logical(M))
  numSigs = length(which(M))
  numPossible = prod(dim(M))
  density = 100*numSigs/numPossible
  nDrugs = dim(M)[1]
  nCells = dim(M)[2]
  nPerDrug = min(apply(M, 1, sum))
  nPerCell = min(apply(M, 2, sum))
  return(sprintf('%d out of %d signatures (density %0.1f%%), covering %d drugs and %d cell lines. \nMin of %d sigs per drug and %d sigs per cell.',
                 numSigs, numPossible, density, nDrugs, nCells, nPerDrug, nPerCell))
}

SummarizeTensor = function(tensor){
  A = !is.na(tensor[,1,])
  SummarizeMatrix(A)
}

CompareTensors = function(T1, T2){
  if(!identical(dim(T1), dim(T2))){
    same = FALSE
  }else{
    d = dimnames(T1)
    same = identical(T1, T2[d[[1]], d[[2]], d[[3]]])
  }
  return(same)
}

LoadCDSigs = function(data=get0('sigs'), debug=FALSE){
  
  baseDir = DataDir('expr/')

  if(debug){
    dataFile = paste0(baseDir, 'drugSigs_smallSample.RData')
    n = 1000
  }else{
    dataFile = paste0(baseDir, 'drugSigs.RData')
    n = NDrugSigs()
  }
  
  if(is.null(data) || nrow(data) != n){
    load(dataFile)
  }else{
    sigs = data
  }
  return(sigs)
}

LoadCDInfo = function(debug=FALSE){
  if(debug){
    load(DataDir('metadata/drugSigInfo_smallSample.RData'))
  }else{
    load(DataDir('metadata/drugSigInfo.RData'))
  }
  return(info)
}

ConstructTensor = function(sigs, info, pThresh=1, specificDose=FALSE, time='all',
                            nCells=NA, nDrugs=NA, cellIds=NA, pertIds=NA, 
                            nPerDrug=1, nPerCell=1, print=TRUE, debug=FALSE,
                            removeDuplicates=TRUE, annot=GetLincsAnnot()){
  
  if(!debug){
    if(nrow(sigs) != NDrugSigs()){warning('Input sigs is unexpected size.')}
    if(nrow(info) != NDrugSigs()){warning('Input info is unexpected size.')}
  }
  
  nGenes = 978
  
  # select data
  info = SelectDataForTensor(info, pThresh=pThresh, specificDose=specificDose, time=time,
                             nCells=nCells, nDrugs=nDrugs, cellIds=cellIds, pertIds=pertIds, 
                             nPerDrug=nPerDrug, nPerCell=nPerCell, print=print, 
                             removeDuplicates=removeDuplicates, annot=annot)
  
  # check that there are at least 2 cell lines
  if(length(unique(info$cell_id)) < 2){
    out = info
    stop('Cannot construct tensor with only one cell type')
  }
  
  # get several statistics on seleted data, mainly at the drug-cell-combination level
  maxP = as.matrix(cast(info[,c('pert_id', 'cell_id', 'pvalue')], pert_id ~ cell_id, max, value='pvalue'))
  maxP[is.infinite(maxP)] = NA
  
  meanP = as.matrix(cast(info[,c('pert_id', 'cell_id', 'pvalue')], pert_id ~ cell_id, mean, value='pvalue'))
  meanP[is.infinite(meanP)] = NA
  
  nAvg = as.matrix(cast(info[,c('pert_id', 'cell_id', 'pvalue')], pert_id ~ cell_id, length, value='pvalue'))
  
  nReps = as.matrix(cast(info[,c('pert_id', 'cell_id', 'replicateCount')], pert_id ~ cell_id, sum, value='replicateCount'))
  
  if(print){print(sprintf('Tensor contains %s', SummarizeMatrix(!is.na(maxP))))}
  
  nSigs = nrow(info)
  
  # filter signatures based on selection
  sigs = sigs[rownames(info),]
  sigs$pert_id = info$pert_id
  sigs$cell_id = info$cell_id
  
  # for each drug/cell pair, compute mean of all available signatures
  if(print){print(sprintf('Computing mean signature per drug-cell pair...'))}
  sigs = data.table(sigs)
  meanSigs  = sigs[, lapply(.SD, mean), by=list(pert_id, cell_id), .SDcols=names(sigs)[1:nGenes]]
  meanSigs = as.data.frame(meanSigs)
  if(print){print(sprintf('...done'))}
  
  # then renormalize to have unit length
  normSigs = as.data.frame(t(apply(meanSigs[,3:(nGenes+2)], 1, function(x){x/Norm2(x)})))
  normSigs$pert_id = meanSigs$pert_id
  normSigs$cell_id = meanSigs$cell_id
  
  # then construct into a tensor
  allCells = unique(normSigs$cell_id)
  allDrugs = unique(normSigs$pert_id)
  
  normSigsC = split(normSigs, f=normSigs$cell_id)
  
  tensor = array(data=NA, dim=c(length(allDrugs), nGenes, length(allCells)),
                 dimnames=list(drugs=allDrugs, genes=GetGeneIdsTensor(), cells=allCells))
  
  for(cell in allCells){
    drugs = normSigsC[[cell]]$pert_id
    stopifnot(!anyDuplicated(drugs))
    A = as.matrix(normSigsC[[cell]][,1:nGenes])
    tensor[drugs,,cell] = A
  }

  # sort drug dimension based on number of signatures available
  if(any(!is.na(pertIds))){
    tensor = tensor[pertIds,,]
  }else{
    tensor = tensor[names(sort(NumSigs(tensor, 'drug'), decreasing = TRUE)),,]
  }
  
  # sort cell type dimension based on number of signatures available
  if(any(!is.na(cellIds))){
    tensor = tensor[,,cellIds]
  }else{
    tensor = tensor[,,names(sort(NumSigs(tensor, 'cell'), decreasing=TRUE))]
  }
  
  return(list(tensor=tensor, maxP=maxP, meanP=meanP, nAvg=nAvg, nReps=nReps, nSigs=nSigs))
}

RestrictToCommonSigs = function(tensorList, print=TRUE){
  # check whether dimnames are identical
  out1 = lapply(tensorList, function(tensor) dimnames(tensor)[[1]]) # drugs
  out2 = lapply(tensorList, function(tensor) dimnames(tensor)[[2]]) # genes
  out3 = lapply(tensorList, function(tensor) dimnames(tensor)[[3]]) # cells
  if(!length(unique(out1))==1){
    drugs = Reduce(intersect, out1)
    tensorList = lapply(tensorList, function(tensor) tensor[drugs,,])
    if(print){print(sprintf('Reducing/reordering to %d drugs', length(drugs)))}
  }
  if(!length(unique(out2))==1){
    genes = Reduce(intersect, out2)
    tensorList = lapply(tensorList, function(tensor) tensor[,genes,])
    if(print){print(sprintf('Reducing/reordering to %d genes', length(genes)))}
  }
  if(!length(unique(out3))==1){
    cells = Reduce(intersect, out3)
    tensorList = lapply(tensorList, function(tensor) tensor[,,cells])
    if(print){print(sprintf('Reducing/reordering to %d cells', length(cells)))}
  }
  
  d = lapply(tensorList, dimnames)
  stopifnot(length(unique(d))==1)
  
  idx = unique(unlist(lapply(tensorList, function(x) which(is.na(x)))))
  return(lapply(tensorList, function(tensor){tensor[idx] = NA; return(tensor)}))
}

WriteTensor2Mat = function(tensor, file){
  if(!exists(file)){
    writeMat(file, T=tensor, pertIds=dimnames(tensor)[[1]], geneIds=dimnames(tensor)[[2]],
             cellIds=dimnames(tensor)[[3]])
  }else{
    warning('Not writing to file as it already exists')
  }
}

SubsetTensor = function(tensor, nDrugs=NA, nCells=NA){
  if(!is.na(nCells)){
    nPerCell = NumSigs(tensor, 'cell')
    cells = names(sort(nPerCell, decreasing=TRUE))[1:nCells]
    tensor = tensor[,,cells]
  }
  if(!is.na(nDrugs)){
    nPerDrug = NumSigs(tensor, 'drug')
    drugs = names(sort(nPerDrug, decreasing = TRUE))[1:nDrugs]
    tensor = tensor[drugs,,]
  }
  return(tensor)
}

NDrugSigs = function(){
  return(201484)
}
