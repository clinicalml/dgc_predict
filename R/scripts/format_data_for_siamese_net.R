# Load drug signatures and corresponding metadata
sigs = LoadCDSigs()
info = LoadCDInfo()

info = SelectDataForTensor(info, pThresh=0.1, specificDose=FALSE,
                           time='all', nPerDrug=3, nPerCell=3)
S = sigs[rownames(info),]


I = info[order(info$cell_id, info$pert_id),]
