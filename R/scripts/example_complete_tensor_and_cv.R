
# Parameters
# Options for method are: mean (1D-Mean), mean2 (2D-Mean), knnd (DNPP), fa_lrtc (Tensor), and others (see matlab/CompleteTensor.m)
method = 'fa_lrtc' 
input_size = 'medium' # small = 5 drugs x 4 genes x 8 cell types, medium = 50 x 978 x 11

# Load data
if(input_size == 'small'){
  load(DataDir('tensors/test/T_test.RData')) 
  tensors = list(meas = T_test)
}else if(input_size == 'medium'){
  load(DataDir('tensors/test/T50_1.RData')) 
  tensors = list(meas = T_meas)
}else{
  stop('input_size not recognized')
}

# Startup matlab
matlab = StartMatlab()

# Run cross-validation <- This isn't working for some reason
#tensors$cv = CrossValidateTensor(matlab, tensor=tensors$meas, methods=method)$tensors

# Complete tensor
tensors$pred = CompleteTensor(matlab, tensors$meas, method)

# ** Be sure to close matlab before you delete the pointer! Otherwise, you may have to reboot!
close(matlab)
rm(matlab)
