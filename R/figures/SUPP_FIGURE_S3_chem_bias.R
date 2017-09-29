
load(ResultsDir('large/entity_specific_accuracy.RData')) 
load(DataDir('metadata/tensor_annot.RData'))

# Compile drug-specific accuracy into a data frame
acc = data.frame(pert_id = annot$pertIds, mean=C$mean$drug, mean2 = C$mean2$drug, 
                 dnpp = C$dnpp$drug, tensor = C$tensor$drug)
rownames(acc) = acc$pert_id

# Load ECFP fingerprints
ECFP = read.table(DataDir('chem_info/ecfp6.tsv'), sep='\t', header=FALSE)
L = GetLincsAnnot()
L = L[!is.na(L$smiles),]
rownames(ECFP) = L$pert_id # ECFP was computed using rcdk directly from L after removing smiles..

# Subset chemical structure and accuracy to tensor pert_ids
E = as.data.frame(ECFP)[acc$pert_id,]
E = na.omit(E)
A = acc[rownames(E),]

# Compute chemical similarity via tanimoto coefficients
n = nrow(E)
TC = array(data=NA, dim=c(n,n))
IDX = apply(E, 1, function(x) as.integer(which(x == 1)))
for(i in 1:n){
  if(i %% 100 == 0){print(sprintf('i=%d', i))}
  for(j in IncreasingSequence(i+1, n)){
    TC[i,j] = JaccardIndex_fromIdx(IDX[[i]], IDX[[j]])
  }
}

# Compute max similarity of each compound with any other compound in the set
TC = SymmetrifyMatrix(TC, diag=NA)
maxTC = apply(TC, 2, max, na.rm=TRUE)
names(maxTC) = rownames(E)
TC = as.data.frame(TC)
rownames(TC) = rownames(E)

# Check this for the 4 compounds highlighted in the paper
idx = which(toupper(annot$pertName) %in% c('HY-11007', 'ABT-751','M-3M3FBS', 'CARBETOCIN'))
perts = annot$pertIds[idx]
L2 = subset(L, pert_id %in% perts)
maxTC[perts]
#
# Output:
# BRD-K91623615 BRD-K09635314 BRD-K97056771 BRD-A32161980 
# 0.2935780     0.3285714     0.3037975     0.4230769 
#
# -> So the good results seen for these compounds cannot simply be explained by having a strong cognate in the dataset.

### Now let's compute overall PCT on subsets restricted by maxTC thresholds
if(~exists('tensors')){
  tensors = LoadTensors(tsize='large', print=TRUE)
}

stopifnot(identical(dimnames(tensors$meas)[[1]], annot$pertIds))
pct = lapply(seq(0.05, 1.0, 0.05), function(thresh){ print(thresh); lapply(tensors$cv, function(tensor) 
  ComputePCT(tensors$meas[names(which(maxTC <= thresh)),,], tensor[names(which(maxTC <= thresh)),,]))})

names(pct) = as.character(seq(0.05, 1.0, 0.05))
mp = melt(pct)
names(mp) = c('PCT', 'method', 'threshold')
mp$threshold = as.numeric(mp$threshold)
mp$method = revalue(mp$method, c('mean'='1D-Mean', 'mean2'='2D-Mean','dnpp'='DNPP', 'tensor'='Tensor'))
ggplot(subset(mp, threshold>=0.3), aes(x=threshold, y=PCT, group=method, color=method)) + 
  geom_line(size=1.4) + theme_bw() + labs(x='maxTC threshold') + theme(text=element_text(size=18)) +
  scale_color_manual(values=unlist(GetMethodColors(longName=TRUE)))
ggsave(PlotDir('chem_structure_bias.pdf'), width=6.5, height=4)
