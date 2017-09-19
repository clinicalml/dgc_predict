load('/Users/rhodos/Desktop/Research/LINCS/data/results/tsize/large/accuracy_per_mode.RData')
load('/Users/rhodos/Desktop/Research/LINCS/data/expr/tensor/tsize/large/tensor_annot.RData')

# Compile drug-specific accuracy into a data frame
acc = data.frame(pert_id = annot$pertIds, mean=C$mean$drug, mean2 = C$mean2$drug, 
                 dnpp = C$dnpp$drug, tensor = C$tensor$drug)

# Load ECFP fingerprints
ECFP = read.table(DataDir('chem_info/ecfp6.tsv'), sep='\t', header=FALSE)
L = GetLincsAnnot()
L = L[!is.na(L$smiles),]
rownames(ECFP) = L$pert_id

# Subset to tensor pert_ids
E = as.data.frame(ECFP)[acc$pert_id,]
E = na.omit(E)

# Compute chemical similarity via tanimoto coefficients
n = nrow(E)
TC = array(data=NA, dim=c(n,n))
IDX = apply(E, 1, function(x) as.integer(which(x == 1)))
for(i in 1:n){
  print(i)
  for(j in IncreasingSequence(i+1, n)){
    TC[i,j] = JaccardIndex_fromIdx(IDX[[i]], IDX[[j]])
  }
}

# Compute max similarity of each compound with any other compound in the set
TC = SymmetrifyMatrix(TC, diag=NA)
maxTC = apply(TC, 2, max, na.rm=TRUE)
names(maxTC) = rownames(E)

# Look at differences in the mean drug-specific accuracy between drugs that have
# a structural cognate, to drugs that don't (where this is defined using various TC thresholds)
methods = c(mean='mean', mean2='mean2', dnpp='dnpp', tensor='tensor')
meanDiff = list()
for(thresh in seq(0.4, 1.0, 0.05)){
  dSim = names(which(maxTC >= thresh))
  acc$hasSim = acc$pert_id %in% dSim
  
  m = melt(acc)
  names(m) = c('pert_id', 'hasSim', 'method', 'PCT')
  print(ggplot(m, aes(x=method, y=PCT, group=interaction(hasSim, method), fill=hasSim)) + 
         geom_boxplot() + ggtitle(sprintf('Thresh = %0.2f', thresh)))
  
  thr = sprintf('%0.2f', thresh)
  meanDiff[[thr]] = lapply(methods, function(mth) mean(subset(acc, hasSim)[[mth]]) - mean(subset(acc, !hasSim)[[mth]]))
}

# Plot simply the differences of means, for each method, across TC thresholds
md = melt(meanDiff)
names(md) = c('PCT_diff', 'method', 'threshold')
ggplot(md, aes(x=threshold, y=PCT_diff, group=method, color=method)) + geom_line()

# Based on this, the largest signal is at a threshold of 0.7. 
thresh = 0.7

# Now let's compute overall PCT on this restricted subset
dUnique = names(which(maxTC < thresh))
#source('/Users/rhodos/Desktop/Research/LINCS/submission/dgc_predict/R/src/DataProc.R')
#tensors = LoadTensors(tsize='large', print=TRUE, loadMergeAndPred = FALSE)

#stopifnot(identical(dimnames(tensors$meas)[[1]] == annot$pertIds))
# pct_all = lapply(tensors$cv, function(tensor) ComputePCT(tensors$))
# ComputePCT(tensors$meas
idx = which(toupper(annot$pertName) %in% c('M-3M3FBS', 'CARBETOCIN', 'ABT-751','GNF-2'))
perts = annot$pertIds[idx]
#L2 = subset(L, pert_id %in% annot$pertIds[idx])
maxTC[perts]
