library(R.matlab)
library(ggplot2)
require(plyr)

writeToFile = TRUE

###### load data
D = list()
for(method in c('mean', 'mean2', 'knnd', 'fa_lrtc')){

    out2 = list()
    for(tensor_size in 5){
      out1 = list()
      for(rep_num in 1:5){
        file = ResultsDir(sprintf('small/obs_density/%s_C_size%d_rep%d.mat',
                       method, tensor_size, rep_num))
        out1[[rep_num]] = readMat(file)$C[,1:6]
      }
      tmp = do.call('rbind', out1)
      dimnames(tmp) = list(replicate=1:25, density=seq(10,60,10))
      out2[[tensor_size]] = as.data.frame(as.table(tmp))
      out2[[tensor_size]]$tensor_size = tensor_size
    }
    df = do.call('rbind', out2)
    df = ChangeColumnName(df, from='Freq', to='corr')

  ### reshape data
  df2 = DataSummary(df[,c('corr', 'density', 'tensor_size')], 
                     varname='corr', groupnames=c('density', 'tensor_size'))
  nDrug = 300
  nGene = 978
  nCell = 15
  n = seq(50, nDrug, length.out = 5)
  df2$num_drugs = as.factor(n[df2$tensor_size])
  
  D[[method]] = df2
}

D2 = ldply(D, data.frame)
df2 = ChangeColumnName(D2, from='.id', to='method')
df2$method = as.factor(df2$method)
df2$method = revalue(df2$method, c('mean'='1D-Mean', 'mean2'='2D-Mean', 'knnd'='DNPP', 'fa_lrtc'='Tensor'))

#if(writeToFile){tiff(file=PlotDir('obs_density.svg'), width=640)}
p = ggplot(df2, aes(x=density, y=corr, color=method, group=method, linetype=method)) +
  geom_errorbar(aes(ymin=corr-sd, ymax=corr+sd), size=1, width=0.1, lty=1, show.legend = FALSE) +
  geom_line(size=2) + 
  xlab("Observation density (%)") +
  ylab("PCT") +
  guides(lty=guide_legend(keywidth=12), color=guide_legend(keywidth=5)) +
  theme_classic() + 
  theme(text=element_text(size=26),
        legend.title = element_blank(),
        legend.position = c(0.75, 0.15)) +
  scale_color_manual(values=unlist(GetMethodColors(longName=TRUE))) +
  scale_linetype_manual(values=unlist(list(`1D-Mean`=2,`2D-Mean`=4,`DNPP`=5,`Tensor`=1))) +
  ylim(c(0.2,0.7))
print(p)
ggsave(PlotDir('obs_density.svg'), height=7, width=8)
#if(writeToFile){dev.off()}

