library(ggplot2)
library(dplyr)
library(plotrix)
library(R.matlab)
library(tidyr)
library(Hmisc)

load(ResultsDir('large/entity_specific_accuracy.RData'))
load(DataDir('metadata/tensor_annot.RData'))

# Get mean gene-gene correlations across all cell types
G = ComputeGeneGeneCor(tensors$meas, cellSpecific = TRUE)
meanG = Reduce("+", G) / length(G)
meanG = apply(array(unlist(G) , c(978,978,length(G))), 1:2, function(x) mean(x, na.rm=TRUE))
diag(meanG) = NA

xm = c('knn','KNN')
color = list(drug = NumSigs(tensors$meas, dim='drug'),
             gene = apply(abs(meanG), 1, function(x) max(x, na.rm=TRUE)),
             cell = NumSigs(tensors$meas, dim='cell'))
leg = list(drug='# cells measured',gene = 'max gene-gene cor', cell = '# drugs measured')
cLim = list(drug = c('yellowgreen','red'),gene = c('red', 'blue'),cell = c('yellowgreen','red'))

palette = list(drug='RdYlGn', gene='RdYlBu', cell='RdYlGn')

size = list(drug=0.5, gene=0.5, cell=1)
alpha = list(drug=0.6, gene=0.9, cell=0.9)
circle = list()

labels = list(drug=annot$pertName, gene='', cell=annot$cellIds)
labels$drug[labels$drug=='HY-11007'] = 'GNF-2'
labels$drug[labels$drug=='M-3M3FBS'] = 'M3'
labelSize = list(drug=8, gene=8, cell=8)
labelRange = list(drug="name %in% c('homoharringtonine', 'terfenadine','dexamethasone', 'JNJ-38877605','M3','GNF-2')", 
                  #y-x > 0.15 | x - y > 0.25 | y > 0.8 | y < 0 | 
                  gene='y > 1',
                  cell="name %in% c('HUH7','NCIH508','SNUC5','HEK293T')")
dodge = list(drug=0, gene=-5, cell=0.1)
xLims = list(drug=range(C[[xm[1]]]$drug, C$tensor$drug, 1, -0.15),
             gene=range(C[[xm[1]]]$gene, C$tensor$gene),
             cell= range(C[[xm[1]]]$cell, C$tensor$cell))
yLims=xLims

for(dim in c('drug', 'gene', 'cell')){
  f = substring(dim, 1, 1)
  xlab = bquote('KNN PCT'[.(f)]) 
  ylab = bquote('Tensor PCT'[.(f)])
  x = rev(C[[xm[1]]][[dim]])
  y = rev(C$tensor[[dim]])
  AUCScatter(x, y, main=sprintf('%s-specific accuracy', capitalize(dim)), 
             color=rev(color[[dim]]), palette=palette[[dim]],
             xlim=xLims[[dim]], ylim=yLims[[dim]],
             labels=rev(labels[[dim]]), labelSize=labelSize[[dim]], 
             labelRange=labelRange[[dim]],dodge=dodge[[dim]],
             xlab=xlab, ylab=ylab, alpha=alpha[[dim]],
             circle=circle[[dim]], cLim = cLim[[dim]], size=size[[dim]],
             legend.position=c(1,0), legend.title=leg[[dim]],
             filename=sprintf('%s_specific_PCT_%s.tiff', dim, xm[2]))
}