
library(rhdf5)
# This differs from the previous version in that I'm using the 6k tensor.


load(DataDir('tensors/tensors6k.RData'))

### Construct features
L = list()

# Average
L$mean$obs  = apply(tensors$meas, c(1,2), mean, na.rm=TRUE)
L$mean$full = apply(tensors$comp, c(1,2), mean)

# Max
L$max$obs  = apply(tensors$meas, c(1,2), max, na.rm=TRUE)
L$max$full = apply(tensors$comp, c(1,2), max)

# Concatenated
#L$allcell$full = UnfoldTensor(tensors$comp, dim=1)

# Cell-specific
for(cell in dimnames(tensors$meas)[[3]]){
  print(cell)
  L[[cell]]$obs = tensors$meas[,,cell]
  L[[cell]]$full = tensors$comp[,,cell]
}

### And do some sanity checks

# For mean and max features, should be (nearly) identical for obs and full when all cell types were measured
idxAllCells = which(NumSigs(tensors$meas, 'drug') == nCells)
stopifnot(max(abs(L$max$obs[idxAllCells,] - L$max$full[idxAllCells,])) < 1e-10)

# Check dimnames (should be identical for everything except 'cat' features)
dn = unlist(lapply(L[setdiff(names(L), c('allcell', 'pca200', 'pca978'))], function(x) lapply(x, function(xx) dimnames(xx))), recursive=FALSE)
stopifnot(all(sapply(dn, function(d) identical(d, dn[[1]]))))

### Write to file
save(L, file=DataDir('expr/drug/tensor_features_for_drug_property_prediction_6ktensor_knn.RData'))
