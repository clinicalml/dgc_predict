

MList = list()

for(cell in c('MDAMB231')){
  for(subset in 'cv'){

    MList = lapply(tensors[[subset]], function(X) cor(X[,,cell], use='pairwise'))
    MList$True = cor(tensors$meas[,,cell], use='pairwise')
    names(MList) = c('1D-Mean', '2D-Mean', 'DNPP', 'Tensor','True')
    MList = rev(MList)
    
    idx_keep = SelectGenesToPlot(MList[[1]], nGenes=200)
    MList = lapply(MList, function(M) M[idx_keep,idx_keep])
    names(dimnames(MList[[1]])) = c('foo','bar')
    p=GMultiHeatmap(MList, clusterRows=TRUE, clusterCols=TRUE, dims=c('gene','cell'),
                    cLim = range(MList, na.rm=TRUE), titles=names(MList), titleSize=24,
                    xlab='', ylab='', colLab=FALSE, rowLab=FALSE, colorMap='PuOr')
    p[[6]] = g_legend(GHeatmap(MList[[1]], dims=c('gene','gene'), legend=TRUE, colorMap='PuOr'))
    tiff(PlotDir(sprintf('gene_correlation_heatmaps_%s_%s.tiff', cell, subset)), width=1200, height=250)
    grid.arrange(grobs=p, layout_matrix=matrix(1:6, ncol=6, byrow=TRUE), widths=c(1,1,1,1,1,0.3))
    dev.off()
  }
}
