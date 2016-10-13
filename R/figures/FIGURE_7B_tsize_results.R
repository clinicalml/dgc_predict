library(R.matlab)

#tsize = list()
PCTf = list()
time = list()
M = list()
methods = c('mean', 'mean2', 'knn', 'fa_lrtc')

tsizes = c('small', 'manydrug', 'manycell', 'large')

for(sz in tsizes){
  print(sz)
#  tsize[[sz]] = list()
  M[[sz]] = list()
  
  #### Load measured tensor
#   if(sz == 'large'){
#     tsize[[sz]]$meas = tensors$meas
#   }else{
#     tsize[[sz]]$meas = LoadTensorMat(DataDir(sprintf('expr/tensor/tsize/%s/%s.mat', sz, sz)))$tensor
#   }
    
  #nD = dim(tsize[[sz]]$meas)[1]
  #nC = dim(tsize[[sz]]$meas)[3]

  #### Load CV results
  file = DataDir(sprintf('results/tsize/%s/%s_tensor_results.mat', sz, sz))
#  if(sz == 'large'){
    library(rhdf5)
    h5ls(file)
    PCTf[[sz]] = h5read(file, 'PCTf')
#     tsize[[sz]]$cv = list(mean=h5read(file, '#refs#/b'),
#                             mean2=h5read(file,'#refs#/c'),
#                             knn=h5read(file,'#refs#/d'),
#                             fa_lrtc=h5read(file,'#refs#/e'))
    time[[sz]] = h5read(file, 'time')
#  }else{
#    out = readMat(file)
#     tsize[[sz]]$cv = list(mean=out$T.imp[[1]][[1]], 
#                       mean2=out$T.imp[[2]][[1]],
#                       knn=out$T.imp[[3]][[1]], 
#                       fa_lrtc=out$T.imp[[4]][[1]])
#    PCTf[[sz]] = out$PCTf
#    time[[sz]] = out$time
#  }
  #tsize[[sz]]$cv$ensemble = (tsize[[sz]]$cv$knn + tsize[[sz]]$cv$fa_lrtc)*(1/2)

  #### Plot PCTf across methods
  stopifnot(nrow(PCTf[[sz]])==10)
  rownames(PCTf[[sz]]) = as.character(1:10)
  listPCTf = list(mean=PCTf[[sz]][,1], mean2=PCTf[[sz]][,2], knn=PCTf[[sz]][,3], fa_lrtc=PCTf[[sz]][,4])
  tiff(PlotDir(sprintf('%s_tensor_results_PCTf.tiff', sz)))
  m = melt(listPCTf)
  names(m) = c('PCT', 'method')
  p = ggplot(m, aes(x=method, y=PCT, fill=method)) +
    stat_summary(fun.y=mean,position=position_dodge(width=0.95),geom="bar") +
    stat_summary(fun.data=mean_sdl, position=position_dodge(0.95),geom="linerange",  size=2) + 
    ggtitle(sprintf('Accuracy on %s tensor', sz)) +
    scale_y_continuous(limits=c(0,1)) + 
    theme_bw() + xlab('') + ylab('PCT per fold') +
    theme(axis.text=element_text(size=20), 
          axis.title=element_text(size=26),
          legend.title = element_blank(),
          legend.text = element_text(size=20), 
          plot.title = element_text(size=26))
  print(p)
  dev.off()
  
  #### Plot time (do it just like PCTf)
  
#   #### Compute individual PCTs
#   M[[sz]] = list()
#   for(method in methods){
#     print(method)
#     M[[sz]][[method]] = matrix(data=NA, nrow=nD, ncol=nC, 
#                                dimnames=list(dimnames(tsize[[sz]]$meas)[[1]], 
#                                              dimnames(tsize[[sz]]$meas)[[3]]))
#     for(i in 1:nD){
#       print(i)
#       for(j in 1:nC){
#         if(!is.na(tsize[[sz]]$meas[i,1,j])){
#           M[[sz]][[method]][i,j] = cor(tsize[[sz]]$cv[[method]][i,,j], tsize[[sz]]$meas[i,,j], use='pairwise')
#         }
#       }
#     }
#   }
#   
#   #### Plot individual PCTs
#   cLim = range(M[[sz]], na.rm = TRUE)
#   for(method in methods){
#     pdf(PlotDir(sprintf('%s_pct_across_%s_tensor.pdf', method, sz)))
#     print(GHeatmap(M[[sz]][[method]], rowLab=FALSE, colLab=FALSE, 
#                    dims=c('drugs', 'cells'), main=toupper(method), cLim=cLim))
#     dev.off()
#   }
}

# out = lapply(PCTf, function(PCT){colnames(PCT) = c('1D-Mean', '2D-Mean', 'KNN', 'Tensor'); return(PCT)})
# pdf(PlotDir('PCTf_all_tensors_all_methods.pdf'), width=12)
# PlotPCTf(out, main='Accuracy per fold across different sized tensors')
# dev.off()

#load(DataDir('results/tsize/tsize_PCTf.RData'))

pdf(PlotDir('PCTf_tsize.pdf'), width=8)
names(PCTf) = c('Small \n(d)','ManyDrugs\n (c)','ManyCells\n (b)','Large\n (a)')
PCTf = lapply(PCTf, function(x){colnames(x) = methods; return(x)})
m = melt(PCTf)
names(m) = c('fold','method', 'PCT','tensor')
m$method = factor(m$method, levels=c('mean', 'mean2', 'knn', 'fa_lrtc'))
m$method = revalue(m$method, c('mean'='1D-Mean', 'mean2'='2D-Mean', 'knn'='KNN', 'fa_lrtc'='Tensor'))
colors = GetMethodColors(longName=TRUE)
p = ggplot(m, aes(x=tensor, y=PCT, fill=method)) +
  stat_summary(fun.y=mean, position=position_dodge(width=0.9),geom="bar", aes(width=0.85)) +
  stat_summary(fun.data=mean_sdl, position=position_dodge(0.9), geom="linerange", size=2) + 
  scale_fill_manual(values=as.character(colors), breaks=names(colors)) +
  #ggtitle('Accuracy evaluated on full tensors') +
  scale_y_continuous(limits=c(0,1)) + 
  theme_classic() + xlab('') + ylab('PCT per fold') +
  theme(axis.text=element_text(size=18), 
        axis.title=element_text(size=26),
        legend.title = element_blank(),
        legend.text = element_text(size=20),
        legend.position = c(0.85, 0.87),
        plot.title = element_text(size=26))
print(p)
dev.off()

# Statistically compare differences between different tensor sizes
names(PCTf) = tsizes
pvals = sapply(methods, function(method) t.test(PCTf$small[,method], PCTf$manycell[,method])$p.value)
if(exists('OUTPUT')){
  OUTPUT$tsize_manycells_v_small_pvals = pvals
  #save(OUTPUT, file=DataDir('results/tsize/OUTPUT.RData'))
}