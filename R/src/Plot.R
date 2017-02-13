

Heatmap = function(C, main=NULL, key=TRUE, file=NULL, margins=c(1,1),
                    colorLim=range(C), colors=c('white', 'blue'), dendrogram='none'){
  nBreaks = 50
 colBreaks = seq(colorLim[1], colorLim[2], length.out = nBreaks+1)
  if(!is.null(file)){
    pdf(file=PlotDir(file=paste0(file, '.pdf')))
  }
  if(all(C==0)){
    warning('all entries of matrix equal 0')
  }
 heatmap.2(C, Rowv=TRUE, Colv=TRUE, trace='none', dendrogram=dendrogram, 
           key=key, main=main, col=colorRampPalette(colors)(nBreaks), 
           breaks=colBreaks, key.title='', margins=margins, labCol='') 
   if(!is.null(file)){
    dev.off()
  }
}

GHeatmap = function(M, dims=names(dimnames(M)), 
                     xlab=dims[1], ylab=dims[2], 
                     cLim=range(M, na.rm=TRUE),
                     clusterRows=FALSE, clusterCols=FALSE, 
                     rowLab=TRUE, colLab=TRUE, xAxisLabSize=12, 
                     labSize=20, titleSize=24,
                     legend=TRUE, main=NULL,
                     low='blue', mid='#444444', high='yellow',
                     colorMap=NA){
  
  # get row and column ordering
  if(clusterRows){
    rowInd = ClusterRows(M)
  }else{
    rowInd = 1:nrow(M)
  }
  
  if(clusterCols){
    colInd = ClusterCols(M)
  }else{
    colInd = 1:ncol(M)
  }
  
  M = t(M[rowInd, colInd])
  
  m = melt(as.matrix(M))
  m[,1] = factor(m[,1], levels=unique(m[,1]), ordered = TRUE)
  m[,2] = factor(m[,2], levels=unique(m[,2]), ordered = TRUE)
  if(!is.null(dims)){
    names(m)[1] = dims[1]
    names(m)[2] = dims[2]
    eval(parse(text=sprintf('p = ggplot(m, aes(x=%s, y=%s, fill=value))', dims[1], dims[2])))
  }else{
    p = ggplot(m, aes(x=Var1, y=Var2, fill=value))
  }
  p = p + geom_tile(aes(fill=value))
  
  if(all(M >= 0, na.rm=TRUE)){
    p = p + scale_fill_gradient(limits=cLim, low = 'White', high = 'NavyBlue')
  }else if(!is.na(colorMap)){
    p = p + scale_fill_gradientn(colours=brewer.pal(7,colorMap))
  }else{
   p = p + scale_fill_gradient2(limits=cLim, low=low, mid=mid, high=high)
  }
  
  p = p + scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    xlab(ylab) + ylab(xlab) +
    theme(axis.title = element_text(size=labSize),
          plot.title = element_text(size=titleSize, hjust=0.5),
          legend.title = element_blank())
  
  if(colLab){
    p = p + theme(axis.text.x = element_text(hjust = 1, angle = 45, size=xAxisLabSize))
  }else{
    p = p + theme(axis.text.x = element_blank(),
                  axis.ticks.x = element_blank())
  }
  
  if(!rowLab){
    p = p + theme(axis.text.y = element_blank(),
                  axis.ticks.y = element_blank())
  }
  
  if(!legend){
    p = p + guides(fill=FALSE)
  }
  
  if(!is.null(main)){
    p = p + ggtitle(main) + theme(plot.title = element_text(hjust = 0.5))
  }
  
  return(p)
}

MultiDens = function(df, main=''){
  df2 = melt(df)
  n = ncol(df)
  p = ggplot(df2, aes(x = value, fill = variable, group=variable)) + 
    geom_density(alpha = 1/n) +
    ggtitle(main) +
    theme(plot.title=element_text(hjust=0.5))
  return(p)
}

SelectGenesToPlot = function(M, nGenes=100){
  idx_var = GetVarianceTopK(M, nGenes)
  idx_med = GetMedianTopK(abs(M), nGenes)
  idx_keep = union(idx_var, idx_med)
  print(sprintf('selecting %d genes', length(idx_keep)))
  return(idx_keep)
}

FancyScatter = function(x, y, xlim=c(0,1), ylim=c(0,1), diag=TRUE, 
                       labels=NULL, labelSize=5, labelRange='2+2==4',
                       color = 'blue', size=1, xlab='', ylab='', circle=NULL, main='',
                       minCount=NULL, maxColor=NULL, cLim = c('red', 'blue'),
                       palette=NULL, print=TRUE, shape=16, alpha=1, legend.position=NULL,
                       legend.title='', filename=NULL, dodge=0.9){
  if(!is.null(minCount)){
    idx = which(color >= minCount)
    x = x[idx]
    y = y[idx]
    color = color[idx]
  }
  
  if(!is.null(maxColor)){
    color[color >= maxColor] = maxColor
  }
  
  df = data.frame(x=x, y=y, color=color, size=size)
  
  p = ggplot(df, aes(x=x, y=y, color=color, size=size, shape=shape, label=labels)) +
    geom_point(alpha=alpha) +
    xlim(xlim) + ylim(ylim) + 
    xlab(xlab) + ylab(ylab) +
    scale_size(guide=FALSE) +
    theme_classic() +
    theme(text = element_text(size=22), 
          legend.background = element_rect(fill='transparent'),
          legend.key.size = unit(1.2,'cm'),
          plot.title = element_text(size=26, hjust=0.5)) + 
    ggtitle(main)

  if(!is.null(circle)){
    p = p + geom_point(aes(x=x[circle], y=y[circle]), shape=21, size=7, color='red') + scale_shape(solid=FALSE)
  }
  
  if(!is.null(cLim) & is.numeric(color)){
    if(!is.null(palette)){
      p = p + scale_color_gradientn(colours=brewer.pal(7,palette), name=legend.title,
                                   guide = guide_colorbar(direction='horizontal', title.position='top'))
    }else{
      p = p + scale_color_gradient(low=cLim[1], high=cLim[2], name=legend.title, 
                                   guide = guide_colorbar(direction='horizontal', title.position='top')) 
    }
  }
  
  if(diag){
    p = p + geom_abline(slope=1, intercept=0) 
  }
  
  if(length(unique(color)) == 1){
    p = p + scale_color_continuous(guide = FALSE)
  }

  if(is.numeric(shape)){
    p = p + scale_shape_identity()
  }
  
  if(!is.null(legend.position)){
    p = p + theme(legend.justification=legend.position, legend.position=legend.position)
  }
  
  if(length(labels) == nrow(df)){
    df$name = labels
    data = eval(parse(text=sprintf('subset(df, %s)', labelRange)))
    p = p + geom_text_repel(data=data, aes(x=x, y=y, group=x, label=name), color='black',
                     size=labelSize, segment.size = 1, segment.color='black', box.padding = unit(1.1, "lines"),
                     point.padding = unit(0.3, "lines"))
    }
  
  if(print){
    print(p)
  }
  
  if(!is.null(filename)){
    tiff(PlotDir(filename), width=600)
    print(p)
    dev.off()
  }
  return(p)
}

CellSpecScatter = function(x, y, method, diag=TRUE,
                            xlim=c(0,1), ylim=c(0,1),
                            xlab='true cell specificity',
                            ylab='cell specificity among predicted profiles',
                            subset){
  if(subset == 'cv'){
    main = 'Only cross-validated predictions'
  }else if(subset == 'pred'){
    main = 'Only true predictions'
  }else if(subset == 'merge'){
    main = 'All predictions'
  }
  
  #main='Preservation of cell specificity'
  df = data.frame(x=x, y=y, method=method)
  
  p = ggplot(df, aes(x=x, y=y, group=method)) +
    geom_point(aes(color=method, shape=method), size=3, alpha=0.6) +
    xlim(xlim) + ylim(ylim) + ggtitle(main) +
    xlab(xlab) + ylab(ylab) +
    theme_classic() +
    theme(legend.text=element_text(size=18),
          legend.position = c(0.16, 0.67),
          legend.title = element_blank(),
          legend.background = element_rect(fill='white', size=0.6, linetype='solid', colour ='DarkGrey'),
          axis.text = element_text(size=20),
          axis.title = element_text(size=20),
          plot.title = element_text(size=20, hjust=0.5),
          plot.margin=unit(c(1,1,1,1),'cm')) +
    scale_color_manual(values=unlist(GetMethodColors(longName = TRUE))) +
    scale_shape_manual(values=c(15, 17, 18, 19)) +
    guides(colour = guide_legend(override.aes = list(size=10)))
  
  if(diag){
    p = p + geom_abline(slope=1, intercept=0, size=1, color='grey', lty=2) 
  }
  
  print(p)
}


PlotCSP = function(cs_true, cs_pred, subset=''){
  nDrug = length(cs_true)
  method = as.factor(rep(c('1D-Mean', '2D-Mean', 'DNPP', 'Tensor'), each=nDrug))
  true = c(rep(cs_true, 4))
  ystr = list(cv='cross-validated', pred='predicted', merge='cv plus predicted')
  pred = c(cs_pred$mean, cs_pred$mean2, cs_pred$dnpp, cs_pred$tensor)
  tiff(PlotDir(sprintf('preservation_of_cell_specificity_%s.tiff', subset)), width=620)
  CellSpecScatter(true, pred, method=method, subset=subset, 
                  xlab = 'cell specificity (measured)',
                  ylab = sprintf('cell specificity (%s)', ystr[[subset]]), 
                  xlim=c(0,1.2), ylim=c(0,1.2))
  dev.off()
}

GetMethodColor = function(method){
  if(method == 'mean'){
    color = 'red'
  } else if(method == 'mean2'){
    color = '#56B4E9'
  }else if(method == 'tensor' || method == 'fa_lrtc'){
    color = '#E69F00'
  }else if(method == 'knnd' || method == 'knn' || method == 'dnpp'){
    color = '#66CC00'
  }else if(method == 'ensemble'){
    color = 'purple'
  }else{
    stop('unexpected method argument')
  }
  return(color)
}

GetMethodColors = function(longName = FALSE){
  colors = list(mean=GetMethodColor('mean'),
              mean2=GetMethodColor('mean2'),
              dnpp = GetMethodColor('dnpp'),
              tensor=GetMethodColor('tensor'))
  if(longName){
    names(colors) = c('1D-Mean', '2D-Mean', 'DNPP', 'Tensor')
  }
  return(colors)
}

PlotPCTf = function(listPCTf, main){
  m = melt(listPCTf)
  names(m) = c('fold', 'method', 'PCT', 'tensor')
  m$method = factor(m$method, levels=c('1D-Mean', '2D-Mean', 'DNPP', 'Tensor', 'Ensemble'))
  colors = GetMethodColors(longName=TRUE)
  p = ggplot(m, aes(x=tensor, y=PCT, fill=method)) +
    stat_summary(fun.y=mean,position=position_dodge(width=0.95),geom="bar") +
    stat_summary(fun.data=mean_cl_normal,position=position_dodge(0.95),geom="linerange",  size=2) + 
    ggtitle(main) +
    scale_fill_manual(values=as.character(colors))+
    scale_y_continuous(limits=c(0,1)) + 
    theme_bw() + xlab('') + ylab('PCT per fold') +
    theme(axis.text=element_text(size=20), 
          axis.title=element_text(size=26),
          legend.title = element_blank(),
          legend.text = element_text(size=20), 
          plot.title = element_text(size=26, hjust=0.5))
  return(p)
}

PlotDataAvailabilityInTensor = function(tensor, main='', xAxisLabSize=6, labSize=28){
  A = MatrixCast(!is.na(tensor[,1,]), type='numeric')
  print(GHeatmap(A, rowLab=FALSE, colLab=TRUE, dims=c('drugs', 'cells'), 
                 main=main, legend=FALSE, xAxisLabSize=xAxisLabSize, labSize=labSize))
}

GMultiHeatmap = function(MList, clusterRows=FALSE, colLab=TRUE, rowLab=TRUE,
                          clusterCols=FALSE, titleSize=24,
                          dims=names(dimnames(M)),xlab=dims[1], ylab=dims[2],
                          cLim = range(MList, na.rm=TRUE), titles=names(MList), 
                          low='blue', mid='#444444', high='yellow',
                          colorMap=NA){
  
  M = MList[[1]]
  
  # get row and column ordering
  if(clusterRows){
    rowInd = ClusterRows(M)
  }else{
    rowInd = 1:nrow(M)
  }
  
  if(clusterCols){
    colInd = ClusterCols(M)
  }else{
    colInd = 1:ncol(M)
  }
  
  MList = lapply(MList, function(M) M[rowInd, colInd])
  pList = lapply(1:length(MList), function(i) 
    GHeatmap(MList[[i]], dims=dims, cLim=cLim, main=titles[i], legend=FALSE, colLab=colLab, rowLab=rowLab,
             xlab=xlab, ylab=ylab, titleSize=titleSize, low=low, mid=mid, high=high, colorMap=colorMap))
  return(pList)
}

g_legend=function(a.gplot){ 
  tmp = ggplot_gtable(ggplot_build(a.gplot)) 
  leg = which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend = tmp$grobs[[leg]] 
  return(legend)} 

PlotGCP = function(GCP, subset='', numSigs, legend=FALSE, title=FALSE, file=NULL){
  if(is.null(file)){
    file = sprintf('preservation_of_ggc_%s.tiff', subset)
  }
  
  if(subset == 'cv'){
    main = 'Only cross-validated predictions'
  }else if(subset == 'pred'){
    main = 'Only true predictions'
  }else if(subset == 'merge'){
    main = 'All predictions'
  }
  
  library(tidyr)
  cellIds = names(GCP[[1]])
  m = list()
  s = list()
  colors = GetMethodColors(longName=TRUE)
  df = data.frame(cell=cellIds, mean=GCP$mean, mean2=GCP$mean2, dnpp=GCP$dnpp, tensor=GCP$tensor)
  df2 = gather(df, 'method', 'cor', 2:5)
  df2$method = revalue(df2$method, c('mean'='1D-Mean', 'mean2'='2D-Mean', 'dnpp'='DNPP', 'tensor'='Tensor'))
  df2 = transform(df2, cell=reorder(as.factor(cell), -numSigs[cell]))
  p = ggplot(df2, aes(x=cell, y=cor, fill = method)) +
    geom_bar(stat='identity', position=position_dodge()) +
    scale_fill_manual(values=as.character(colors), breaks=names(colors)) +
    scale_y_continuous(limits=c(0,1)) + 
    theme_bw() + 
    theme(axis.text.x=element_text(size=14, angle=30, hjust=0.8), 
          axis.text.y=element_text(size=15),
          axis.title=element_text(size=20)) +
    xlab('') + ylab('GCP')
  
  if(legend){
    p = p + theme(legend.title = element_blank(),legend.text = element_text(size=16))
  }else{
    p = p + theme(legend.position='none')
  }
  
  if(title){
    p = p + ggtitle(main) + theme(plot.title = element_text(size=30, hjust=0.5))
  }
    
  tiff(PlotDir(file), width=2000, height=300)
  print(p)
  dev.off()
  m[[subset]] = apply(as.matrix(df[,2:4]), 2, mean)
  s[[subset]] = apply(as.matrix(df[,2:4]), 2, sd)
}


SubsetROC = function(roc.in, n_short){
  if(length(roc.in) > 2*n_short){
    n_long = length(roc.in$tpr)
    idx = SubsetROCIdx(n_long, n_short)
    roc.out = list()
    roc.out$tpr = roc.in$tpr[idx]
    roc.out$fpr = roc.in$fpr[idx]
  }else{
    roc.out = roc.in
  }
  return(roc.out)
}

# choose indices between 1 and n_long so that I have approximately n_short
# indices evenly spaced and including the first and last points
SubsetROCIdx = function(n_long, n_short){
  n_skip = floor(n_long / n_short)
  idx = seq(1, n_long, n_skip)
  if( idx[length(idx)] != n_long){
    idx = c(idx, n_long)
  }
  return(idx)
}

