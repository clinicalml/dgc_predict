
### Determine which method has better performance per signature
PCT = laply(tensors$cv, function(tensor) ComputePCTPerSig(tensors$meas, tensor, format='list')$R)
dimnames(PCT) = list(method = names(tensors$cv), 
                     drug = dimnames(tensors$meas)[[1]],
                     cell = dimnames(tensors$meas)[[3]])

M = apply(PCT, c(2,3), function(x){a = names(x)[which.max(x)]; return(ifelse(length(a)==0, NA, a))})

colors = GetMethodColors(longName=FALSE)
p = GHeatmap(M, rowLab=FALSE, colLab=TRUE, dims=c('drugs', 'cells'),
             ylab = 'cell lines', main='', 
             legend=FALSE, xAxisLabSize=5.5, labSize=28) + 
  scale_fill_manual(values=colors)
print(p)
ggsave(file=PlotDir('best_method_per_drug_cell.eps'), width=7, height=7)
