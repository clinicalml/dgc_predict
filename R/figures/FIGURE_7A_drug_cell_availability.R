
tsize = 'large'

### Load measured tensor
tensors = list(meas=LoadTensorMat(DataDir(sprintf('tensors/%s.mat', tsize)))$tensor)

### Load cross-validated tensors
file = ResultsDir(sprintf('%s/%s_tensor_results.mat', tsize, tsize))
tensors$cv = list(mean=h5read(file, '#refs#/b'), mean2=h5read(file,'#refs#/c'),
                  dnpp=h5read(file,'#refs#/d'), tensor=h5read(file,'#refs#/e'))
#tensors$cv$ensemble = 0.5 * (tensors$cv$dnpp + tensors$cv$tensor)
tensors$cv = lapply(tensors$cv, function(tensor){dimnames(tensor) = dimnames(tensors$meas); return(tensor)})

### Determine which method has better performance per signature
PCT = laply(tensors$cv, function(tensor) ComputePCTPerSig(tensors$meas, tensor, format='list')$R)
dimnames(PCT) = list(method = names(tensors$cv), 
                     drug = dimnames(tensors$meas)[[1]],
                     cell = dimnames(tensors$meas)[[3]])

M = apply(PCT, c(2,3), function(x){a = names(x)[which.max(x)]; return(ifelse(length(a)==0, NA, a))})

colors = GetMethodColors(longName=FALSE)
p = GHeatmap(M, rowLab=FALSE, colLab=TRUE, dims=c('drugs', 'cells'),
             ylab = 'cell lines', main='', 
             legend=FALSE, xAxisLabSize=7, labSize=28) + 
  scale_fill_manual(values=colors)
print(p)
ggsave(file=PlotDir('best_method_per_drug_cell.png'), width=9, height=9)



