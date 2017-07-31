
L = GetLincsAnnot()
load(DataDir('metadata/tensor_annot_6k.RData'))

A = L[L$pert_id %in% annot[[1]],]
targets_list = Str2Vec(na.omit(A$targets), split=':')
targets = unique(unlist(targets_list))
tCount = sapply(targets, function(tar){length(which(sapply(
  targets_list, function(tlist) tar %in% tlist)))})
myTargets = names(tCount[tCount >= 20])

# Construct matrix
pertsWithTargets = A$pert_id[!is.na(A$targets)]
nPerts = length(pertsWithTargets)

targetMatrix = array(data=0, dim=c(nPerts, length(myTargets)),
                     dimnames=list(perts=pertsWithTargets, targets=myTargets))

tList = Str2Vec(A$targets, split=':')
names(tList) = A$pert_id


for(target in myTargets){
  pertIds = names(which(sapply(tList, function(myT) target %in% myT)))
  targetMatrix[pertIds,target] = 1
}

### Sanity check
a = apply(targetMatrix, 2, sum)
b = tCount[tCount >= 20]
stopifnot(all(a==b))

write.table(targetMatrix, file=DataDir('targets/top_targets_in_6ktensor.csv'), 
            sep=',', col.names=TRUE, row.names=TRUE)
