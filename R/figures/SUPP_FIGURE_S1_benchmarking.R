library(R.matlab)
library(reshape2)
library(ggplot2)
library(Rmisc)

# LOAD DATA 
results_list = list()
for(model in c('asmatrix', 'constrained', 'si_lrtc', 'ha_lrtc', 'fa_lrtc','tmac','mean','mean2','knnd')){
  print(model)
  out = readMat(ResultsDir(sprintf('small/benchmarking/%s_full_small_tensor.mat', model)))
  results_list[[model]] = list(PCTf=out$PCTf, time=out$time)
}
results = melt(results_list)[,3:5]
names(results) = c('value', 'type', 'model')

results$model = revalue(as.factor(results$model), c('mean'='1D-Mean', 'mean2'='2D-Mean', 'knnd'='DNPP',
                                                    'fa_lrtc'='FaLRTC','si_lrtc'='SiLRTC','ha_lrtc'='HaLRTC',
                                                    'asmatrix'='AsMatrix','constrained'='Constrained','tmac'='Tmac'))
text_size = 16
res1 = subset(results, type=='PCTf')
res1 = transform(res1, model=reorder(model, value))
p1 = ggplot(res1) + theme_bw() + theme(text = element_text(size=text_size), plot.title=element_text(hjust=0.5)) +
  geom_boxplot(aes(x=model, y=value, group=model)) + ggtitle('Accuracy') + xlab('') + ylab('PCT per fold') +
  theme(axis.text.x = element_text(hjust = 1, angle = 45), text=element_text(size=20))

res2 = subset(results, type=='time')
res2$model = reorder(res2$model,as.numeric(res1$model),order=TRUE)
p2 = ggplot(res2) + theme_bw() + theme(text = element_text(size=text_size), plot.title=element_text(hjust=0.5)) +
  geom_boxplot(aes(x=model, y=value, group=model)) + ggtitle('Runtime') + xlab('') + ylab('seconds per fold') +
  scale_y_log10() + theme(axis.text.x = element_text(hjust = 1, angle = 45), text=element_text(size=20))

pdf(PlotDir('benchmarking_results.pdf'), width=12)
multiplot(p1, p2, cols=2)
dev.off()


       
