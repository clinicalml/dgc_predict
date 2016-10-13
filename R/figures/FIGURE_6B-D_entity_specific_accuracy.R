library(ggplot2)
library(dplyr)
library(plotrix)
library(R.matlab)
library(tidyr)
library(Hmisc)

# Get accuracy per mode
#source('R/analyze_tensor_results/compute_accuracy_per_mode.R')
load(DataDir('results/tsize/large/accuracy_per_mode.RData'))

# Get tensor annotations
load(DataDir('/expr/tensor/tsize/large/tensor_annot.RData'))

# Get gene gene correlation matrix
# (the next line assumes you have already computed the cell-specific correlations)
# G = ComputeGeneGeneCor(tensors$meas, cellSpecific = TRUE)
#meanG = Reduce("+", G) / length(G)
#meanG = apply(array(unlist(G) , c(978,978,length(G))), 1:2, function(x) mean(x, na.rm=TRUE))
#diag(meanG) = NA

xm = c('knn','KNN')
#xm = c('mean','1D-Mean')

color = list(drug = NumSigs(tensors$meas, dim='drug'),
             gene = apply(abs(meanG), 1, function(x) max(x, na.rm=TRUE)),
             cell = NumSigs(tensors$meas, dim='cell'))
leg = list(drug='# cells measured',gene = 'max gene-gene cor', cell = '# drugs measured')
cLim = list(drug = c('yellowgreen','red'),gene = c('red', 'blue'),cell = c('yellowgreen','red'))

palette = list(drug='RdYlGn', gene='RdYlBu', cell='RdYlGn')#list(drug='YlOrBr', gene='YlGnBu', cell='YlOrBr')

size = list(drug=0.5, gene=0.5, cell=50) # hm doesn't seem to have any effect
alpha = list(drug=0.6, gene=0.9, cell=0.9)
circle = list()

labels = list(drug=annot$pertName, gene='', cell=annot$cellIds)
labels$drug[labels$drug=='HY-11007'] = 'GNF-2'
labels$drug[labels$drug=='M-3M3FBS'] = 'M3'
labelSize = list(drug=8, gene=8, cell=8)
labelRange = list(drug="name %in% c('homoharringtonine', 'terfenadine','dexamethasone', 'JNJ-38877605','M3','GNF-2')", 
                  #y-x > 0.15 | x - y > 0.25 | y > 0.8 | y < 0 | 
                  gene='y > 1', #0.75 | y < 0.2', 
                  cell="name %in% c('HUH7','NCIH508','SNUC5','HEK293T')")
                  #cell='y <=1')
#'y < 0.42 | y > 0.64 | abs(y-x) > 0.1' 
dodge = list(drug=0, gene=-5, cell=0.1)
xLims = list(drug=range(C[[xm[1]]]$drug, C$tensor$drug, 1, -0.15),#C$mean$drug, C$mean2$drug, 
             gene=range(C[[xm[1]]]$gene, C$tensor$gene), #C$mean$gene, C$mean2$gene, 
             cell= range(C[[xm[1]]]$cell, C$tensor$cell)) #C$mean$cell, C$mean2$cell, 
#xLims = list(drug=c(-0.15,1), gene=c(-0.15,1), cell=c(-0.15,1))
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

# # Compute correlation between PCTd and cell-specificity
# load(DataDir('results/tsize/large/cell_specificity.RData'))
# lapply(list('1D-Mean'='mean','2D-Mean'='mean2','KNN'='knn','Tensor'='tensor'),
#        function(method) cor(C[[method]]$drug, cs$true))