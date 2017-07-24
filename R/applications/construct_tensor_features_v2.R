
library(rhdf5)
# This differs from the previous version in that I'm not mapping to pubchem Ids

### Load measured tensor
tensors = list(meas=LoadTensorMat(DataDir('tensors/large.mat'))$tensor)
names(dimnames(tensors$meas)) = c('drug','gene','cell')

### Load completed tensors
comp = list()
for(method in c('dnpp', 'tensor')){
  print(method)
  tcomp = h5read(ResultsDir(sprintf('large/hdf5/%s_final_hdf5.mat', method)),'T')
  dimnames(tcomp) = dimnames(tensors$meas)
  comp[[method]] = NormSigs(tcomp)
}

### Compute ensemble mean
tensors$comp = 0.5*(comp$dnpp + comp$tensor)

### Map entrez ids to gene symbols because who likes to look at numeric id's all day long
dimnames(tensors$meas)[[2]] = MapEntrezToSymbol(dimnames(tensors$meas)[[2]], lm=TRUE)
dimnames(tensors$comp)[[2]] = MapEntrezToSymbol(dimnames(tensors$comp)[[2]], lm=TRUE)

### Also subset to top ten measured cell types
tensors$meas = tensors$meas[,,1:10]
tensors$comp = tensors$comp[,,1:10]

# and then remove drugs that have no signatures among these ten
idx = which(NumSigs(tensors$meas, 'drug')==0)
tensors$meas = tensors$meas[-idx,,]
tensors$comp = tensors$comp[-idx,,]

# # Write pertIds and corresponding pubchemIds to file
# pertIds = annot$pertIds[match(pubchemIds, annot$pubchemIds)]
# drugAnnot = data.frame(pertID=pertIds, pubchemID=pubchemIds)
# write.table(drugAnnot, file=DataDir('metadata/drugAnnot_for_HaoChih.txt'),
#                        row.names=FALSE, col.names=TRUE, quote=FALSE)

### Construct features
X = list()

# Average
X$mean_obs  = apply(tensors$meas, c(1,2), mean, na.rm=TRUE)
X$mean_full = apply(tensors$comp, c(1,2), mean)

# Max
X$max_obs  = apply(tensors$meas, c(1,2), max, na.rm=TRUE)
X$max_full = apply(tensors$comp, c(1,2), max)

# Concatenated
X$cat_full = UnfoldTensor(tensors$comp, dim=1)

# Cell-specific
for(cell in dimnames(tensors$meas)[[3]]){
  print(cell)
  X[[paste(cell, 'meas', sep='_')]] = tensors$meas[,,cell]
  X[[paste(cell, 'full', sep='_')]] = tensors$comp[,,cell]
}

### And do some sanity checks

# for mean and max features, should be (nearly) identical for obs and full when all cell types were measured
idxAllCells = which(NumSigs(tensors$meas, 'drug') == 10)
stopifnot(max(abs(X$max_obs[idxAllCells,] - X$max_full[idxAllCells,])) < 1e-10)

### Write to file
save(X, file=DataDir('side_effect/tensor_features_for_side_effect_prediction.RData'))
