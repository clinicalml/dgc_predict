load(DataDir('cell_viability/CTRP.RData'))
set.seed(123)

### Shuffle things just in case there was a reason for the input ordering
stopifnot(rownames(GR50) == rownames(GRmax))
reorder = sample(rownames(GR50), size=nrow(GR50), replace=FALSE)
GR50 = GR50[reorder,]
GRmax = GRmax[reorder,]
GR50_meas = GR50_meas[reorder,]
GRmax_meas = GRmax_meas[reorder,]
stopifnot(rownames(GR50) == rownames(GRmax))
stopifnot(rownames(GR50_meas) == rownames(GR50))

### Load large tensors, measured and predicted
#tensors = LoadTensors(tsize='large', loadMergeAndPred=TRUE)
#tensors = list(meas=tensors$meas, comp=tensors$merge)

### Hold out data from (the same) 22 drugs (10%) as test data for both GR50 and GRmax
n = nrow(GRmax)
num_obs_per_drug = apply(GRmax_meas, 1, function(x) length(which(!is.na(x))))
num_obs_per_cell = apply(GRmax_meas, 2,  function(x) length(which(!is.na(x))))
#idx_7_or_8_obs = which(num_obs_per_drug >= 7)
#idx_holdout = sample(idx_7_or_8_obs, size=22, replace = FALSE)
#idx_keep = setdiff(1:n, idx_holdout)
# Try holding out 22 random drugs, and seeing how many this corresponds to in the actual data per cell line
idx_holdout = sample(1:n, size=22, replace=FALSE)
Ymeas_50 = GR50_meas[idx_holdout,]
Ymeas_M = GRmax_meas[idx_holdout,]
num_obs_per_cell_50 = apply(Ymeas_50, 2, function(x) length(which(!is.na(x))))
num_obs_per_cell_M = apply(Ymeas_M, 2, function(x) length(which(!is.na(x))))

idx_keep = setdiff(1:n, idx_holdout)

Ymeas_50 = GR50_meas[idx_keep,]
Ymeas_M = GRmax_meas[idx_keep,]
stopifnot(rownames(Ymeas_50) == rownames(Ymeas_M))
num_obs_per_cell_50 = apply(Ymeas_50, 2, function(x) length(which(!is.na(x))))
num_obs_per_cell_M = apply(Ymeas_M, 2, function(x) length(which(!is.na(x))))

Y_50 = GR50[idx_keep,]
Y_M = GRmax[idx_keep,]
Ytest_50 = GR50[idx_holdout,]
Ytest_M = GRmax[idx_holdout,]


### Define 6 different data inputs:


# Baselines (based on measured signatures only):
# 1) Cell-matched signatures, only measured ***
# 2) Averaged signatures across most common cell lines, only measured

# Our inputs:
# 3) Cell-matched signatures, measured + imputed ***
# 4) Averaged signatures across most common cell lines, measured + imputed
# 5) All cell lines
# 6) Drug embedding

### Setup cell-matched signatures, only measured

### Use 3-fold CV to train and test LASSO models to predict cell viability <- will extend to random forests, SVM, and maybe DNNs.

# For each of the 8 cell lines


