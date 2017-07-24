
# Load CTRP data and append pert_ids
CTRP  = read.table(DataDir('cell_viability/CTRP_GRmetrics.tsv'), sep = '\t', header=TRUE)
drugs = unique(CTRP$Drug)
drugs = setdiff(drugs, drugs[grepl(':', drugs)])
ctrpMap = ReadXLS(DataDir('cell_viability/CTRPv2.0._INFORMER_SET.xlsx')) # mapping from CTRP drug names to BRD ids
ctrpAnnot = data.frame(name=drugs, pert_id = ctrpMap$broad_cpd_id[match(drugs, ctrpMap$cpd_name)])
CTRP = subset(CTRP, Drug %in% drugs)
CTRP$PertID = ctrpAnnot$pert_id[match(CTRP$Drug, ctrpAnnot$name)]
CTRP$Cell.line[CTRP$Cell.line == 'U266B1'] = 'U266'

# Intersect with tensor info
load(DataDir('metadata/tensor_annot.RData'))
commonPerts = intersect(annot$pertIds, ctrpAnnot$pert_id) # 220 perturbagens in common
commonCells = intersect(annot$cellIds, unique(CTRP$Cell.line)) # 7 cell lines in common

# Subset to intersection
CTRP = subset(CTRP, PertID %in% commonPerts)
CTRP = subset(CTRP, Cell.line %in% commonCells)

# Extract two matrices from CTRP 
library(reshape)
GR50 = cast(CTRP, PertID ~ Cell.line, value='GR50')
perts = GR50$PertID
cells = colnames(GR50)[-1]
GR50 = as.matrix(GR50[,-1])
dimnames(GR50) = list(perts=perts, cells=cells)
GR50[is.infinite(GR50)] = NA

GRmax = cast(CTRP, PertID ~ Cell.line, value='GRmax')
rownames(GRmax) = GRmax$PertID
GRmax = GRmax[,-1]

stopifnot(identical(rownames(GR50), rownames(GRmax)))
stopifnot(identical(colnames(GR50), colnames(GRmax)))

# Count available data that maps to MEASURED L1000 signatures in tensor
# tensor = LoadTensorMat(DataDir('tensors/large.mat'))$tensor
# tSlice = tensor[,1,]
S = tSlice[rownames(GRmax), colnames(GRmax)]
idxMissing = which(is.na(S), arr.ind=TRUE)

GR50_meas = GR50
GR50_meas[idxMissing] = NA
GRmax_meas = GRmax
GRmax_meas[idxMissing] = NA


# Count available data per cell line
overall_counts = count(CTRP$Cell.line)
ct_meas = apply(as.matrix(GR50_meas), 2, function(x) length(which(!is.na(x))))
counts = overall_counts
counts$ct_meas = ct_meas
counts = ChangeColumnName(counts, from=c('x', 'freq'), to=c('cell', 'ct_total'))
nSigs = NumSigs(tensors$meas, 'cell')
counts$sig_total = nSigs[counts$cell]

# Plot
mcounts = melt(counts)
ggplot(mcounts, aes(x=cell, y=value, fill=variable)) + 
  geom_bar(stat='identity', position='dodge') +
  ylab('Number of Examples') +
  ggtitle('Data availability in L1000 Drug Tensor / CTRP Intersection')
ggsave(PlotDir('data_counts.pdf'))

GHeatmap(log2(GR50), rowLab=FALSE, low='red', mid='white', high='green', xlab='', ylab='', main='log2(GR50)')
ggsave(PlotDir('GR50_heatmap.pdf'))

G = as.matrix(GRmax)
dimnames(G) = dimnames(GR50)
GHeatmap(G, rowLab=FALSE, low='red', mid='white', high='green', xlab='', ylab='', main='GRmax')
ggsave(PlotDir('GRmax_heatmap.pdf'))

# Let's look at pairwise correlations between these different cell types
library(GGally)
ggpairs(GRmax)
ggsave(PlotDir('GRmax_pairwise_cors.pdf'))

GR50_mat = data.frame(as.matrix(GR50))
ggpairs(GR50_mat)
ggsave(PlotDir('GR50_pairwise_cors.pdf'))

# Some correlations are quite high (0.88) and all are positive. Many are weak.
# Interesting patterns in the scatter plots of GRmax, indicating (IMO) important
# differences between cell types for various drugs.
save(GR50, GRmax, GR50_meas, GRmax_meas, CTRP, file=DataDir('cell_viability/CTRP.RData'))
