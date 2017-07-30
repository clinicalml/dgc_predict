
# I will just identify the top 6 ATC codes associated with drugs in the large tensor
L = GetLincsAnnot()
load('~/projects/side_info/data/annot/tensor_annot.RData')
A = L[L$pert_id %in% annot$pertIds,]

atc = Str2Vec(na.omit(A$atc_code), split='[|]')
atc_letters = lapply(atc, function(ATC) sapply(ATC, function(y) substr(y, 1, 1)))
all_letters = sort(unique(unlist(atc_letters)))

n = list()
for(letter in all_letters){
  n[[letter]] = length(which(sapply(atc_letters, function(ATC) letter %in% ATC)))
}

myATC = c('D','C','L','N','J','A')


# Construct matrix
pertsWithATC= A$pert_id[!is.na(A$atc_code)]
nPerts = length(pertsWithATC)

ATCMatrix = array(data=0, dim=c(nPerts, length(myATC)),dimnames=list(perts=pertsWithATC, atc=myATC))

atcList = Str2Vec(A$atc_code, split='[|]')
names(atcList) = A$pert_id
just_letters = lapply(atcList, function(ATC) sapply(ATC, function(y) substr(y, 1, 1)))

for(letter in myATC){
  pertIds = names(which(sapply(just_letters, function(myA) letter %in% myA)))
  ATCMatrix[pertIds,letter] = 1
}

write.table(ATCMatrix, file=DataDir('atc/top_6_atc_codes_in_large_tensor.csv'), 
            sep=',', col.names=TRUE, row.names=TRUE)
