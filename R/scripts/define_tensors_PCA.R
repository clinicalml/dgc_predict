library(tensor)

## TODO: 
# TEST UNCENTER
# TEST COMPUTE ERROR RATE
# RETHINK HOW I'M COMPRESSING AND EVALUATING

Uncenter = function(T_est, bias, std){
  return(aperm(laply(1:length(std), function(i) T_est[,i,]*std[i] + bias[i]), c(2, 1, 3)))
}

# Set parameters
maxPC = 978
pcs = c(seq(10, 100, 10), seq(200, maxPC, 100), maxPC)
print = TRUE

# Load drug signatures and corresponding metadata
sigs = LoadCDSigs()
info = LoadCDInfo()

# Subset metadata and get corresponding signatures
info_select = SelectDataForTensor(info, pThresh=0.1, specificDose=FALSE, removeDuplicates=TRUE,
                               time='all', nPerDrug=3, nPerCell=3, print=TRUE)
S = sigs[rownames(info_select),]


# Center data
bias = apply(S, 2, mean)
std = apply(S, 2, sd)
S = t(laply(1:length(std), function(i) (S[,i] - bias[i])/std[i]))

# Compute (UNREGULARIZED) SVD <- This might be okay in our case since we have n >> p (n = roughly 25*p)
list[d, U, V] = svd(S, nu=maxPC, nv=maxPC)
U = as.data.frame(U)
names(U) = paste0('PC', 1:maxPC)

# Append drug and cell info to factor matrix
U$pert_id = info_select$pert_id
U$cell_id = info_select$cell_id

# For each drug/cell pair, compute mean of all available signatures
if(print){print(sprintf('Computing mean signature per drug-cell pair...'))}
U = data.table(U)
meanSigs  = as.data.frame(U[, lapply(.SD, mean), by=list(pert_id, cell_id), .SDcols=names(U)[1:maxPC]])

# # then renormalize to have unit length
# normSigs = as.data.frame(t(apply(meanSigs[,3:(maxPC+2)], 1, function(x){x/Norm2(x)})))
# normSigs$pert_id = meanSigs$pert_id
# normSigs$cell_id = meanSigs$cell_id

#(actually this version is not normalized)
normSigs = data.frame(meanSigs[,3:(maxPC+2)], pert_id=meanSigs$pert_id, cell_id=meanSigs$cell_id)

# then construct into a tensor
allCells = unique(normSigs$cell_id)
allDrugs = unique(normSigs$pert_id)

normSigsC = split(normSigs, f=normSigs$cell_id)

tensor = array(data=NA, dim=c(length(allDrugs), maxPC, length(allCells)),
               dimnames=list(drugs=allDrugs, genes=paste0('PC', 1:maxPC), cells=allCells))

for(cell in allCells){
  drugs = normSigsC[[cell]]$pert_id
  stopifnot(!anyDuplicated(drugs))
  A = as.matrix(normSigsC[[cell]][,1:maxPC])
  tensor[drugs,,cell] = A
}

# sort drug and cell type dimensions based on number of signatures available
tensor = tensor[names(sort(NumSigs(tensor, 'drug'), decreasing = TRUE)),,]
tensor = tensor[,,names(sort(NumSigs(tensor, 'cell'), decreasing=TRUE))]

# define 4 tensors
tensors = list()
tensors$large = tensor
tensors$manycell = SubsetTensor(tensor, nDrugs=300, nCells=50)
tensors$manydrug = SubsetTensor(tensor, nDrugs=2000, nCells=15)
tensors$small = tensor[dimnames(tensors$manycell)[[1]],,dimnames(tensors$manydrug)[[3]]] 
stopifnot(all(sapply(tensors, is.numeric)))

# Startup matlab
matlab = StartMatlab()

# For each dimension in range, run cross-validation and project back into gene space
T_meas = LoadTensorMat(DataDir('tensors/small.mat'))$tensor
PCT = list()
PCTs = list()
PCT_compressed = list()
L2_compressed = list()
for(nPC in 3:30){
  print(sprintf('******* Starting calculations for %d principal components ********', nPC))
  npc = as.character(nPC)

  print('Running cross-validation...')
  out = CrossValidateTensor(matlab, tensor=tensors$small[,1:nPC,])
  tensors$pca[[npc]] = out$tensors
  PCT_compressed[[npc]] = out$PCTf
  
  L2_compressed[[npc]] = lapply(out$tensors, function(tensor) Tensor2Vec(tensors$small, tensor)

  print('Projecting back to the full space...')
  DV = diag(d[1:nPC]) %*% t(V[,1:nPC])
  tensors$full[[npc]] = lapply(tensors$pca[[npc]], function(tensor){print('one tensor..'); aperm(tensor(tensor, DV, 2, 1), c(3,1,2))})

  print('Computing overall PCT...')
  PCT[[npc]] = lapply(tensors$full[[npc]], function(T_est) ComputePCT(T_meas, Uncenter(aperm(T_est, c(2,1,3)), bias=bias, std=std)))

  print('Computing PCT per signature...')
  PCTs[[npc]] = lapply(tensors$full[[npc]], function(T_est) ComputePCTPerSig(T_meas, Uncenter(aperm(T_est, c(2,1,3)), bias=bias, std=std)))
}

# Shut down matlab
close(matlab)
rm(matlab)

# Plot overall PCT as a function of the number of principal components
mPCT = melt(PCT)
qplot(data=mPCT, x=as.numeric(L1), y=value, group=L2, color=L2, xlab='# of PCA dimensions', ylab='PCT (overall accuracy)') + 
  geom_line(size=1.2) + theme(axis.title=element_text(size=18))
ggsave(file=PlotDir('overall_PCT_vs_PCA_dimension_V2.pdf'))

# Plot PCT in compressed space vs. number of principal components
m = melt(PCT_compressed)
names(m) = c('fold','method','value','nPC')
m$nPC = as.numeric(m$nPC)
m = subset(m, nPC <= 30)
m$method = revalue(m$method, c('mean'='1D-Mean', 'mean2'='2D-Mean', 'fa_lrtc'='Tensor', 'knnd'='DNPP'))
#M = DataSummary(mPCTo, varname='value', groupnames=c('nPC', 'method'))
qplot(data=m, x=as.numeric(nPC), y=value, group=method, color=method, xlab='# of PCA dimensions', ylab='PCT per fold') + 
  theme(axis.title=element_text(size=18)) + 
  geom_point(size=0.1, alpha=0.2) + geom_smooth()
ggsave(file=PlotDir('PCT_compressed_vs_PCA_dimension_30.pdf'))

# Plot PCT per signature as a function of the number of principal components
mPCTs = subset(melt(PCTs), variable=='R')
names(mPCTs) = c('variable', 'value', 'method', 'nPC')
M = DataSummary(mPCTs, varname='value', groupnames=c('nPC', 'method'))
M$method = revalue(M$method, c('mean'='1D-Mean', 'mean2'='2D-Mean', 'fa_lrtc'='Tensor', 'knnd'='DNPP'))
qplot(data=M, x=as.numeric(nPC), y=value, group=method, color=method, xlab='# of PCA dimensions', ylab='PCT per signature') + 
  geom_line(size=1.2) + theme(axis.title=element_text(size=18)) + 
  geom_errorbar(aes(ymin=value-sd, ymax=value+sd), size=0.5, position=position_dodge(10))
ggsave(file=PlotDir('PCT_per_signature_vs_PCA_dimension_V2.pdf'))

tsrs = list(meas=tensors$small, full=tensors$full, pca=tensors$pca)
save(PCT, PCTs, tsrs, file=PlotDir('PCA_results.RData'))
