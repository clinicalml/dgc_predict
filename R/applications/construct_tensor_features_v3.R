
library(rhdf5)
# This differs from the previous version in that I'm adding the 'cv' tensor.
# This is so that I can compare the use of measured vs. predicted profiles on
# the same set.

tsize = 'large'

### Load measured tensor
tensors = list(meas=LoadTensorMat(DataDir(sprintf('tensors/%s.mat', tsize)))$tensor)
names(dimnames(tensors$meas)) = c('drug','gene','cell')

### Load completed tensors
comp = list()
for(method in c('dnpp', 'tensor')){
  print(method)
  tcomp = h5read(ResultsDir(sprintf('%s/hdf5/%s_final_hdf5.mat', tsize, method)),'T')
  dimnames(tcomp) = dimnames(tensors$meas)
  comp[[method]] = NormSigs(tcomp)
}

### Load cross-validated tensors
file = ResultsDir(sprintf('%s/%s_tensor_results.mat', tsize, tsize))
tcv = list(dnpp=h5read(file,'#refs#/d'), tensor=h5read(file,'#refs#/e'))
tcv = lapply(tcv, function(tensor){dimnames(tensor) = dimnames(tensors$meas); return(tensor)})

### Compute ensemble mean
tensors$comp = 0.5*(comp$dnpp + comp$tensor)
tensors$cv = 0.5*(tcv$dnpp + tcv$tensor)

### Map entrez ids to gene symbols because who likes to look at numeric id's all day long
dimnames(tensors$meas)[[2]] = MapEntrezToSymbol(dimnames(tensors$meas)[[2]], lm=TRUE)
dimnames(tensors$comp)[[2]] = MapEntrezToSymbol(dimnames(tensors$comp)[[2]], lm=TRUE)
dimnames(tensors$cv)[[2]] = MapEntrezToSymbol(dimnames(tensors$cv)[[2]], lm=TRUE)

### Also subset to top ten measured cell types
tensors$meas = tensors$meas[,,1:10]
tensors$comp = tensors$comp[,,1:10]
tensors$cv = tensors$cv[,,1:10]

# and then remove drugs that have no signatures among these ten
idx = which(NumSigs(tensors$meas, 'drug')==0)
tensors$meas = tensors$meas[-idx,,]
tensors$comp = tensors$comp[-idx,,]
tensors$cv = tensors$cv[-idx,,]

### Construct features
X = list()

# Average
X$mean$obs  = apply(tensors$meas, c(1,2), mean, na.rm=TRUE)
X$mean$cv  = apply(tensors$cv, c(1,2), mean, na.rm=TRUE)
X$mean$full = apply(tensors$comp, c(1,2), mean)

# Max
X$max$obs  = apply(tensors$meas, c(1,2), max, na.rm=TRUE)
X$max$cv  = apply(tensors$cv, c(1,2), max, na.rm=TRUE)
X$max$full = apply(tensors$comp, c(1,2), max)

# Concatenated
X$cat$full = UnfoldTensor(tensors$comp, dim=1)

# Cell-specific
for(cell in dimnames(tensors$meas)[[3]]){
  print(cell)
  X[[cell]]$obs = tensors$meas[,,cell]
  X[[cell]]$cv = tensors$cv[,,cell]
  X[[cell]]$full = tensors$comp[,,cell]
}

### And do some sanity checks

# For mean and max features, should be (nearly) identical for obs and full when all cell types were measured
idxAllCells = which(NumSigs(tensors$meas, 'drug') == 10)
stopifnot(max(abs(X$max_obs[idxAllCells,] - X$max_full[idxAllCells,])) < 1e-10)

# Check obs pattern: should be (nearly the) same between meas and cv
for(nm in names(X)){
  if(names(X[[nm]]) == c('obs','cv','full')){
    idx1 = which(is.na(X[[nm]]$obs))
    idx2 = which(is.na(X[[nm]]$cv))
    print(sprintf('%s: %d sigs differ', nm, length(union(setdiff(idx1, idx2), setdiff(idx2, idx1)))/978))
  }else{
    print(sprintf('%s does not contain all versions of features', nm))
  }
}

# Check dimnames (should be identical for everything except 'cat' features)
dn = unlist(lapply(X[setdiff(names(X), 'cat')], function(x) lapply(x, function(xx) dimnames(xx))), recursive=FALSE)
stopifnot(all(sapply(dn, function(d) identical(d, dn[[1]]))))

# Check correlation between observed and cv for matched vs. unmatched indices, two examples
cor(X$PC3$obs[1,], X$PC3$cv[1,]) #matched drug
cor(X$PC3$obs[2,], X$PC3$cv[1,]) #un-matched drug
cor(X$VCAP$obs[,2], X$VCAP$cv[,2], use = 'pairwise') #matched gene
cor(X$VCAP$obs[,100], X$VCAP$cv[,3], use = 'pairwise') #un-matched gene

### Write to file
save(X, file=DataDir('expr/drug/tensor_features_for_drug_property_prediction.RData'))
