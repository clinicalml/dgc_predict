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
# save(nAvg, meanP, maxP, nReps, file=DataDir('metadata/large_tensor_stats.RData'))
load(DataDir('metadata/tensor_stats.RData'))

#### Load cmap data
data = readMat(DataDir('tensors/cmap.mat'))
CMAP = aperm(data$C, c(2, 3, 1))
dimnames(CMAP) = list(drug=unlist(data$dim2C), gene=unlist(data$dim3C), cell=unlist(data$dim1C))
cmapNames = dimnames(CMAP)[[1]]

#### Load mapping between LINCS and CMAP drugs
# Load mapping
load(DataDir('metadata/pairwiseMap.RData'))
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

#### Finally, compute the correlations
C = lapply(tensors$cv, function(tensor) ComputePCTPerSig(CMAP, SubsetTensorBy(LINCS, tensor), format='df'))
names(C) = c('1D-Mean', '2D-Mean', 'DNPP', 'Tensor')
C$True = ComputePCTPerSig(CMAP, LINCS, format='df')
sig = C$True$adjP < 0.05 & C$True$R > 0
for(method in names(C)){
  C[[method]]$method = method
  C[[method]]$sig = sig
}

allC = Reduce(rbind, C)
allC$group = 'all'
allC_sig = subset(allC, sig)
allC_sig$group = 'significant'

data = rbind(allC, allC_sig)

data$method = factor(data$method, levels=c('True', '1D-Mean', '2D-Mean', 'DNPP', 'Tensor'), ordered = TRUE)

p = ggplot(data, aes(x=method, y=R, group=interaction(group,method), fill=group)) +
  geom_boxplot() +
  theme_bw() +
  xlab('Method') +
  ylab('Correlation per signature') +
  theme(text = element_text(size=28), legend.title = element_blank(), 
        legend.text=element_text(size=28), legend.position = c(0.85, 0.1)) +
  scale_fill_brewer(palette='YlGnBu') +
  ylim(c(-0.6,1))

tiff(PlotDir('CMAP_L1000_comparisons.tiff'), width=600, height=500)
print(p)
dev.off()

diff = lapply(C, function(c) t.test(C$True[sig,'R'], c[sig,'R'], paired=TRUE)$estimate)
pval = lapply(C, function(c) t.test(C$True[sig,'R'], c[sig,'R'], paired=TRUE)$p.value)
# Turns out that DNPP is the only one you can't distinguish from the correlations with measured signatures
