
library(rhdf5)

tsize = 'large'

### Load measured tensor
tensors = list(meas=LoadTensorMat(DataDir(sprintf('tensors/%s.mat', tsize)))$tensor)
names(dimnames(tensors$meas)) = c('drug','gene','cell')

### Load completed tensor
tensors$comp = h5read(ResultsDir(sprintf('%s/hdf5/dnpp_final_hdf5.mat', tsize)),'T')
dimnames(tensors$comp) = dimnames(tensors$meas)

### Load cross-validated tensor
file = ResultsDir(sprintf('%s/%s_tensor_results.mat', tsize, tsize))
tensors$cv = h5read(file,'#refs#/d')
dimnames(tensors$cv) = dimnames(tensors$meas)

### Map entrez ids to gene symbols 
dimnames(tensors$meas)[[2]] = MapEntrezToSymbol(dimnames(tensors$meas)[[2]], lm=TRUE)
dimnames(tensors$comp)[[2]] = MapEntrezToSymbol(dimnames(tensors$comp)[[2]], lm=TRUE)
dimnames(tensors$cv)[[2]] = MapEntrezToSymbol(dimnames(tensors$cv)[[2]], lm=TRUE)

### Also subset to top ten measured cell types
nCells = 10
tensors$meas = tensors$meas[,,1:nCells]
tensors$comp = tensors$comp[,,1:nCells]
tensors$cv = tensors$cv[,,1:nCells]

# and then remove drugs that have no signatures among these ten
idx = which(NumSigs(tensors$meas, 'drug')==0)
tensors$meas = tensors$meas[-idx,,]
tensors$comp = tensors$comp[-idx,,]
tensors$cv = tensors$cv[-idx,,]

### Construct features
L = list()

# Average
L$mean$obs  = apply(tensors$meas, c(1,2), mean, na.rm=TRUE)
L$mean$cv  = apply(tensors$cv, c(1,2), mean, na.rm=TRUE)
L$mean$full = apply(tensors$comp, c(1,2), mean)

# Max
L$max$obs  = apply(tensors$meas, c(1,2), max, na.rm=TRUE)
L$max$cv  = apply(tensors$cv, c(1,2), max, na.rm=TRUE)
L$max$full = apply(tensors$comp, c(1,2), max)

# Concatenated
L$allcell$full = UnfoldTensor(tensors$comp, dim=1)

# Concatenated with reduced dimensionality
out = prcomp(L$allcell$full, retx=TRUE, center=TRUE, scale.=FALSE)
cs = cumsum(out$sdev^2/sum(out$sdev^2))
plot(cs)
L$pca200$full = out$x[,1:200] # The first 200 accounts for 64.0% of the variance
L$pca978$full = out$x[,1:978] # The first 978 accounts for 90.7% of the variance

# Cell-specific
for(cell in dimnames(tensors$meas)[[3]]){
  print(cell)
  L[[cell]]$obs = tensors$meas[,,cell]
  L[[cell]]$cv = tensors$cv[,,cell]
  L[[cell]]$full = tensors$comp[,,cell]
}

### And do some sanity checks

# For mean and max features, should be (nearly) identical for obs and full when all cell types were measured
idxAllCells = which(NumSigs(tensors$meas, 'drug') == nCells)
stopifnot(max(abs(L$max$obs[idxAllCells,] - L$max$full[idxAllCells,])) < 1e-10)

# Check obs pattern: should be (nearly the) same between meas and cv
for(nm in names(L)){
  if(names(L[[nm]]) == c('obs','cv','full')){
    idx1 = which(is.na(L[[nm]]$obs))
    idx2 = which(is.na(L[[nm]]$cv))
    print(sprintf('%s: %d sigs differ', nm, length(union(setdiff(idx1, idx2), setdiff(idx2, idx1)))/978))
  }else{
    print(sprintf('%s does not contain all versions of features', nm))
  }
}

# Check dimnames (should be identical for everything except 'cat' features)
dn = unlist(lapply(L[setdiff(names(L), c('allcell', 'pca200', 'pca978'))], function(x) lapply(x, function(xx) dimnames(xx))), recursive=FALSE)
stopifnot(all(sapply(dn, function(d) identical(d, dn[[1]]))))

# Check correlation between observed and cv for matched vs. unmatched indices, two examples
cor(L$PC3$obs[1,], L$PC3$cv[1,]) #matched drug
cor(L$PC3$obs[2,], L$PC3$cv[1,]) #un-matched drug
cor(L$VCAP$obs[,2], L$VCAP$cv[,2], use = 'pairwise') #matched gene
cor(L$VCAP$obs[,100], L$VCAP$cv[,3], use = 'pairwise') #un-matched gene

### Write to file
save(L, file=DataDir(sprintf('expr/tensor_features_for_drug_property_prediction_%dcells_knn.RData', nCells)))
