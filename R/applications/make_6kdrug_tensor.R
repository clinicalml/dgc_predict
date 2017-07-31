# Load drug signatures and corresponding metadata
sigs = LoadCDSigs()
info = LoadCDInfo()

# Construct tensor with as many drugs as possible, filtering by pThresh=0.1
out = ConstructTensor(sigs=sigs, info=info, pThresh=0.1, specificDose=FALSE, time='all',
                          print=TRUE, nPerDrug=1, nPerCell=1, removeDuplicates=TRUE, nCells=10)

dimnames(out$tensor)
NumSigs(out$tensor, 'cell')
NumSigs(out$tensor, 'drug')

tensors = list()
tensors$meas = out$tensor

matlab = StartMatlab()
tensors$comp = CompleteTensor(matlab, tensors$meas, method='knnd')

dimnames(tensors$meas)[[2]] = MapEntrezToSymbol(dimnames(tensors$meas)[[2]], lm=TRUE)
dimnames(tensors$comp)[[2]] = MapEntrezToSymbol(dimnames(tensors$comp)[[2]], lm=TRUE)

save(tensors, file=DataDir('tensors/tensors6k.RData'))
annot = dimnames(tensors$meas)
save(annot, file=DataDir('metadata/tensor_annot_6k.RData'))
