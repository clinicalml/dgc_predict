

library(easyGgplot2)
library(RColorBrewer)

MList = list()

for(cell in c('MDAMB231')){
  for(subset in 'cv'){
#     MList = list(True=G$meas[[cell]], Tensor=G[[subset]]$tensor[[cell]], KNN=G[[subset]]$knn[[cell]], 
#                 '2D-Mean'=G[[subset]]$mean2[[cell]],'1D-Mean'=G[[subset]]$mean[[cell]])

    MList = lapply(tensors[[subset]], function(X) cor(X[,,cell], use='pairwise'))
    MList$True = cor(tensors$meas[,,cell], use='pairwise')
    names(MList) = c('1D-Mean', '2D-Mean', 'KNN', 'Tensor','True')
    MList = rev(MList)
    
    #drugKeep = which(!is.na(tensors$meas[,1,cell]))
    #idx_keep = SelectGenesToPlot(t(tensors$meas[drugKeep,,cell]), nGenes=200)
    
    idx_keep = SelectGenesToPlot(MList[[1]], nGenes=200)
    MList = lapply(MList, function(M) M[idx_keep,idx_keep])
    names(dimnames(MList[[1]])) = c('foo','bar')
    p=GMultiHeatmap(MList, clusterRows=TRUE, clusterCols=TRUE, dims=c('gene','cell'),
                    cLim = range(MList, na.rm=TRUE), titles=names(MList), titleSize=24,
                    xlab='', ylab='', colLab=FALSE, rowLab=FALSE)
    p[[6]] = g_legend(GHeatmap(MList[[1]], dims=c('gene','gene'), legend=TRUE))
    tiff(PlotDir(sprintf('gene_correlation_heatmaps_%s_%s.tiff', cell, subset)), width=1200, height=250)
    grid.arrange(grobs=p, layout_matrix=matrix(1:6, ncol=6, byrow=TRUE), widths=c(1,1,1,1,1,0.3))
    dev.off()
  }
}