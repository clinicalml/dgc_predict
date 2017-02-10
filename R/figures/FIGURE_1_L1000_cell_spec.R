load(ResultsDir('large/cell_specificity.RData'))

# print distribtuion across 2130 drugs
tiff(PlotDir('cell_specificity_L1000.tiff'), height=300, width=480)
data = data.frame(cs=cs$true)
p = ggplot(data, aes(x=cs, color='black')) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  stat_density(alpha=.2, fill="#FF6600", bw=0.04) +
  theme(text=element_text(size=24),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank()) +
  xlab('Cell specificity') +
  ylab('Density in L1000')
print(p)
dev.off()

# choose 4 examples

annot = GetLincsAnnot()
rownames(data) = dimnames(tensors$meas)[[1]]
data$nSig = NumSigs(tensors$meas, 'drug')
data$name = annot$name[match(rownames(data), annot$pert_id)]
d2 = subset(data, nSig >= 5)

s = seq(0.38,1.05,length.out=4)
idx = lapply(1:4, function(i) which(abs(d2$cs-s[i]) < 0.02))
low = which(d2$name == 'HOMOHARRINGTONINE') #0.38
mid1 = which(d2$name == 'TERFENADINE') # cs=0.61
mid2 = which(d2$name == 'DEXAMETHASONE') #cs=0.83
high = which(d2$name == 'JNJ-38877605') #cs=1.05

size=47
for(example in c(low, mid1, mid2, high)){
  pert_id = rownames(d2)[example]
  name = d2$name[example]
  M = GetDrugSlice(tensors$meas, pert_id)
  idx_keep = SelectGenesToPlot(M, nGenes=40)
  p = GHeatmap(M[idx_keep,], dims=c('gene', 'cell'), clusterRows=TRUE, clusterCols=TRUE,
               main='', xAxisLabSize=60, labSize=60, rowLab=FALSE, legend=FALSE,
               colLab=FALSE)
  if(example != low){
    p = p + xlab('') + ylab('')
  }
  p = p + theme(text=element_text(size=size),
                axis.ticks = element_blank())
                
  tiff(PlotDir(sprintf('%s.tiff', name)))
  print(p)
  dev.off()
}
