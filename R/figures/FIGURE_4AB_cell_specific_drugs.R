
library(gridExtra)


GetDrugSlice <- function(tensor, drug){
  return(t(na.omit(t(tensor[drug,,]))))
}

annot = GetLincsAnnot()
idx = which(annot$name %in% c('ABT-751','M-3M3FBS','HY-11007'))

for(i in idx){
  name = annot$name[i]
  pert = annot$pert_id[i]
  if(pert %in% dimnames(tensors$meas)[[1]]){
    print(name)
    MList = lapply(tensors$cv, function(tensor) GetDrugSlice(tensor, pert))
    MList$meas = GetDrugSlice(tensors$meas, pert)
    idx_keep = SelectGenesToPlot(MList$meas, nGenes=78)
    MList = lapply(MList, function(M) M[idx_keep,])
    MList = rev(MList)
    names(MList) = c('True', 'Tensor', 'KNN', '2D-Mean', '1D-Mean')
    tiff(PlotDir(sprintf('%s.tiff', name)), width=1200, height=350)
    pList = GMultiHeatmap(MList, clusterRows=TRUE, clusterCols=TRUE, dims=c('gene','cell'), 
                          colLab=TRUE, rowLab=FALSE, xlab='', ylab='', titleSize=24)
    pList[[6]] = g_legend(GHeatmap(MList[[1]], dims=c('gene','cell'), legend=TRUE))
    grid.arrange(grobs=pList, left = textGrob('gene', gp=gpar(fontsize=28,font=8), rot=90, hjust=0.2),
                 layout_matrix = matrix(1:6, ncol=6, byrow=TRUE), widths=c(1,1,1,1,1,0.2))
    dev.off()
  }
}
