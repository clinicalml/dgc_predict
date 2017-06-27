# Load drug signatures and corresponding metadata
sigs = LoadCDSigs()
info = LoadCDInfo()

# Construct large tensor
out = ConstructTensor(sigs=sigs, info=info, pThresh=0.1, specificDose=FALSE, time='all',
                      print=TRUE, nPerDrug = 3, nPerCell=3, removeDuplicates=TRUE)

# define tensor and corresponding dimnames
manydrug = SubsetTensor(out$tensor, nDrugs=2000, nCells=15)
dimnames = dimnames(manydrug)

# save manydrug tensor to python pickle file
library(rPython)
wd = '~/Desktop/Research/TC2_Lincs/data/ManyDrug/tensor/'
setwd(wd)
python.assign('wd', wd)
python.exec( 'import pickle' )
python.exec( 'import os' )
outName = paste0('manydrug')
manyDrug_noname = unname(manydrug, force = TRUE)
python.assign(outName, manyDrug_noname)
command = paste0("pickle.dump(", outName, ", open('", outName, ".p', 'wb'))")
python.exec(command)
write.table(dimnames[[1]], file='pert_ids.txt', sep='\n', quote=FALSE, col.names=FALSE, row.names=FALSE)
write.table(dimnames[[2]], file='gene_ids.txt', sep='\n', quote=FALSE, col.names=FALSE, row.names=FALSE)
write.table(dimnames[[3]], file='cell_ids.txt', sep='\n', quote=FALSE, col.names=FALSE, row.names=FALSE)
