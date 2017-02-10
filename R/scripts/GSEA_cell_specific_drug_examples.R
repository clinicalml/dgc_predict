GetDrugSlice = function(tensor, drug){
  return(t(na.omit(t(tensor[drug,,]))))
}

distfunc = function(x) as.dist((1-cor(x))/2)

set.seed(42)
annot = GetLincsAnnot()
gsea = list()
results = list()
minGeneSetSize = 5


# --------------------
# M3

pert = annot$pert_id[which(annot$name == 'M-3M3FBS')]
M = GetDrugSlice(tensors$meas, pert)
cells1 = c('A549','AGS','RKO','MCF7')
cells2 = setdiff(colnames(M), union(cells1, c('SKB', 'VCAP')))
stopifnot(length(cells1)==4)
stopifnot(length(cells2)==9)
sig_sm = apply(M[,cells1], 1, mean)
sig_lg = apply(M[,cells2], 1, mean)

gsea$M$sm = RunGsea(profile=sig_sm, dz_genes_up = names(sig_sm), 
                         minGeneSetSize=minGeneSetSize, doGSOA=FALSE, doGSEA=TRUE)
gsea$M$lg = RunGsea(profile=sig_lg, dz_genes_up = names(sig_sm), 
                         minGeneSetSize=minGeneSetSize, doGSOA=FALSE, doGSEA=TRUE)
results$M3 = lapply(gsea$M, function(x) GetGSEAResultsV2(x, merge=TRUE, analysis='GSEA', colname='Adjusted.Pvalue', thresh=0.05))


# --------------------
# HY-11007 (AKA GNF-2)

pert = annot$pert_id[which(annot$name == 'HY-11007')]
M = GetDrugSlice(tensors$meas, pert)
cells2 = setdiff(colnames(M), 'ASC')
stopifnot(length(cells2)==4)
sig_sm = M[,'ASC']
sig_lg = apply(M[,cells2], 1, mean)

gsea$H$sm = RunGsea(profile=sig_sm, dz_genes_up = names(sig_sm), 
                         minGeneSetSize=minGeneSetSize, doGSOA=FALSE, doGSEA=TRUE)
gsea$H$lg = RunGsea(profile=sig_lg, dz_genes_up = names(sig_sm), 
                         minGeneSetSize=minGeneSetSize, doGSOA=FALSE, doGSEA=TRUE)
results$HY = lapply(gsea$H, function(x) GetGSEAResultsV2(x, merge=TRUE, analysis='GSEA', colname='Adjusted.Pvalue', thresh=0.05))


# --------------------
# Carbetocin
set.seed(84)
pert = annot$pert_id[which(annot$name == 'CARBETOCIN')]
M = GetDrugSlice(tensors$meas, pert)


for(cell in colnames(M)){
  print(cell)
  gsea$C[[cell]] = RunGsea(profile=M[,cell], dz_genes_up = rownames(M), 
                           minGeneSetSize=minGeneSetSize, doGSOA=FALSE, doGSEA=TRUE)
}
results$CT = lapply(gsea$C2, function(x) GetGSEAResultsV2(x, merge=TRUE, analysis='GSEA', colname='Adjusted.Pvalue', thresh=0.05))

save(results, file=ResultsDir('large/gsea_cell_specific_drugs.RData'))
