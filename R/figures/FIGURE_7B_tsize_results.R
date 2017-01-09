library(R.matlab)
library(rhdf5)

PCTf = list()
M = list()
methods = c('mean', 'mean2', 'knn', 'fa_lrtc')

tsizes = c('small', 'manydrug', 'manycell', 'large')

for(sz in tsizes){
  print(sz)
  M[[sz]] = list()
  
  #### Load CV results
  file = ResultsDir(sprintf('%s/%s_tensor_results.mat', sz, sz))
  h5ls(file)
  PCTf[[sz]] = h5read(file, 'PCTf')

}
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
