
# Compute cell specificity of each drug in each tensor for Figure 4C
# for(subset in c('cv', 'pred','merge')){
#  list[OUTPUT$CSp_cor[[subset]], OUTPUT$CSp_slope[[subset]], cs[[subset]]] = 
#     ComputeCSpPres(tensors$meas, tensors[[subset]], plot=TRUE, subset=subset, cs_true=cs$true, debug=FALSE)
# }
#save(cs, file=ResultsDir('large/cell_specificity.RData'))

load(ResultsDir('large/cell_specificity.RData'))

# remove results for 5 drugs that had signatures that were not cross-validated (just due to uneven fold split)
df = as.data.frame(cs$pred)
idx = which(df$mean > 1e-15) # the cell-specificity in 1D-Mean 'true predictions' should be 0
cs_pred = lapply(cs$pred, function(x) return(x[-idx]))
cs_true = cs$true[-idx]

PlotCSP(cs$true, cs$cv, subset='cv')
PlotCSP(cs_true, cs_pred, subset='pred')
PlotCSP(cs$true, cs$merge, subset='merge')

# Compute mean squared error as reported in paper
# MSE = list()
# methods=c('mean'='mean','mean2'='mean2','dnpp'='dnpp','tensor'='tensor')
# for(subset in c('cv','pred','merge')){
#   MSE[[subset]] = lapply(methods, function(method) mean((cs$true - cs[[subset]][[method]])^2))
# }

