library(biomaRt)

#### Rerun ConstructTensor in order to get meanP and nReps
# sigs = LoadCDSigs()
# info = LoadCDInfo()
# out = ConstructTensor(sigs=sigs, info=info, pThresh=0.1, specificDose=FALSE, time='all',
#                       nCells=72, nDrugs=10000, print=TRUE, nPerDrug=3, nPerCell=3, removeDuplicates=TRUE)
# meanP = out$meanP
# maxP = out$maxP
# nReps = out$nReps
# nAvg = out$nAvg
# save(nAvg, meanP, maxP, nReps, file=DataDir('expr/tensor/tsize/large/tensor_stats.RData'))
load(DataDir('expr/tensor/tsize/large/tensor_stats.RData'))

#### Load cmap data (This is the data as prepared by the Probabilistic Connnectivity Mapping paper; 718 drugs measured in all three cell lines)
data = readMat('/Users/rhodos/Desktop/Dropbox/Rachel/Thesis/drug_repurposing/data/cmap/tensor_from_ProbCMap.mat')
CMAP = aperm(data$C, c(2, 3, 1))
dimnames(CMAP) = list(drug=unlist(data$dim2C), gene=unlist(data$dim3C), cell=unlist(data$dim1C))
cmapNames = dimnames(CMAP)[[1]]

#### Load mapping between LINCS and CMAP drugs
# Load mapping
load('/Users/rhodos/Desktop/Dropbox/Rachel/Thesis/drug_repurposing/data/lincs/mapping/pairwiseMap.RData')
drugMap = pairwiseMap
drugMap$num.match = sapply(Str2Vec(drugMap$matching.fields, split='[,]'), length)
pertIds = dimnames(tensors$meas)[[1]]

# Restrict the mapping to compounds that are in the tensor and the CMAP subset
drugMap = drugMap[which(drugMap$lincs.id %in% pertIds),] 
drugMap = drugMap[which(drugMap$cmap.name %in% cmapNames),]

# Make sure the mapping is 1:1
drugMap = ddply(drugMap, 'lincs.id', summarise, 
                cmap.name=cmap.name[which.max(num.match)], 
                num.match=max(num.match))
drugMap = ddply(drugMap, 'cmap.name', summarise, 
                lincs.id=lincs.id[which.max(num.match)], 
                num.match=max(num.match))

#### Restrict both tensors to the drugs that map and the appropriate cell lines
CMAP = CMAP[cmapNames %in% drugMap$cmap.name,,c('PC3', 'MCF7','HL60')]
LINCS = tensors$meas[pertIds %in% drugMap$lincs.id,,c('PC3', 'MCF7','HL60')]

# change the drugnames of CMAP to pertIds, and reorder according to LINCS ordering
dimnames(CMAP)[[1]] = drugMap$lincs.id[match(dimnames(CMAP)[[1]], drugMap$cmap.name)]
CMAP = CMAP[dimnames(LINCS)[[1]],,]

# Some drugs don't have data in any of the three cell types in LINCS, so restrict further
idxZero = which(NumSigs(LINCS,'drug')==0)
LINCS = LINCS[-idxZero,,]
CMAP = CMAP[-idxZero,,]

### Get mapping between ensemble gene ids and entrez ids
ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org",
                  path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl")
geneIds = dimnames(LINCS)[[2]]
geneMap = getBM(attributes=c('ensembl_gene_id', 'entrezgene'), 
                filters='entrezgene', values=geneIds, mart=ensembl)
geneMap = geneMap[geneMap$ensembl_gene_id %in% dimnames(CMAP)[[2]],]
geneMap = RemoveDuplicateRowsMulti(geneMap, colnames(geneMap))

# change the gene names of CMAP to entrez ids, and reorder/subset LINCS accordingly
CMAP = CMAP[,geneMap$ensembl_gene_id,]
dimnames(CMAP)[[2]] = as.character(geneMap$entrezgene)
LINCS = LINCS[,as.character(geneMap$entrezgene),]

SubsetTensorBy <- function(LINCS, tensor){
  return(tensor[dimnames(LINCS)[[1]], dimnames(LINCS)[[2]], dimnames(LINCS)[[3]]])
}

#### Finally, compute the correlations
C = lapply(tensors$cv, function(tensor) ComputePCTPerSig(CMAP, SubsetTensorBy(LINCS, tensor), format='df'))
names(C) = c('1D-Mean', '2D-Mean', 'KNN', 'Tensor')
C$True = ComputePCTPerSig(CMAP, LINCS, format='df')
sig = C$True$adjP < 0.05 & C$True$R > 0
for(method in names(C)){
  C[[method]]$method = method
  C[[method]]$sig = sig
}

OUTPUT$CMAP_COMPARISON$ALL = lapply(C, function(x) c(mean=mean(x$R), sd=sd(x$R)))
OUTPUT$CMAP_COMPARISON$SIG = lapply(C, function(x) c(mean=mean(x$R[sig]), sd=sd(x$R[sig])))
OUTPUT$CMAP_COMPARISON$num_higher$all = sapply(C, function(c) length(which(C$True$R - c$R < 0)))
OUTPUT$CMAP_COMPARISON$num_higher$sig = sapply(C, function(c) length(which(C$True$R[sig] - c$R[sig] < 0))) 

allC = Reduce(rbind, C)
allC$group = 'all'
allC_sig = subset(allC, sig)
allC_sig$group = 'significant'

data = rbind(allC, allC_sig)

data$method = factor(data$method, levels=c('True', '1D-Mean', '2D-Mean', 'KNN', 'Tensor'), ordered = TRUE)

p = ggplot(data, aes(x=method, y=R, group=interaction(group,method), fill=group)) +
  geom_boxplot() +
  theme_bw() +
  xlab('Method') +
  ylab('Correlation per signature') +
  #ggtitle('Correlation between CMAP (measured) and \nL1000 (measured or predicted) signatures') +
  theme(text = element_text(size=28), legend.title = element_blank(), 
        legend.text=element_text(size=28), legend.position = c(0.85, 0.1)) +
  scale_fill_brewer(palette='YlGnBu') +
  ylim(c(-0.6,1))

tiff(PlotDir('CMAP_L1000_comparisons.tiff'), width=600, height=500)
print(p)
dev.off()

diff = lapply(C, function(c) t.test(C$True[sig,'R'], c[sig,'R'], paired=TRUE)$estimate)
pval = lapply(C, function(c) t.test(C$True[sig,'R'], c[sig,'R'], paired=TRUE)$p.value)
# Turns out that KNN is the only one you can't distinguish from the correlations with measured signatures


meanp = as.vector(meanP[dimnames(LINCS)[[1]], dimnames(LINCS)[[3]]])
AUCScatter(na.omit(as.numeric(C$True$R))[sig], na.omit(as.numeric(C$KNN$R))[sig], xlim=c(-1,1), ylim=c(-1,1), 
           color=na.omit(-log10(meanp))[sig], labels='')

### Now perform similar analysis to see if intrinsic (i.e. CV) benchmark distinguishes between these methods
C = lapply(tensors$cv, function(tensor) ComputePCTPerSig(LINCS, SubsetTensorBy(LINCS, tensor), format='df'))
names(C) = c('1D-Mean', '2D-Mean', 'KNN', 'Tensor')
for(method in names(C)){
  C[[method]]$method = method
  C[[method]]$sig = sig ## I am keeping the same sig as defined earlier
}

allC = Reduce(rbind, C)
allC$group = 'all'
allC_sig = subset(allC, sig)
allC_sig$group = 'significant'

data = rbind(allC, allC_sig)

data$method = factor(data$method, levels=c('1D-Mean', '2D-Mean', 'KNN', 'Tensor'), ordered = TRUE)

p = ggplot(data, aes(x=method, y=R, group=interaction(group,method), fill=group)) +
  geom_boxplot() +
  theme_bw() +
  xlab('Method') +
  ylab('Correlation per signature') +
  theme(text = element_text(size=28), legend.title = element_blank(), 
        legend.text=element_text(size=28), legend.position = c(0.85, 0.1)) +
  scale_fill_brewer(palette='YlGnBu') +
  ylim(c(-0.6,1))

tiff(PlotDir('L1000_internal_comparisons_on_CMAP_subset.tiff'), width=600, height=500)
print(p)
dev.off()
