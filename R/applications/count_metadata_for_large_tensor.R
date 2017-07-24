load('/Users/rhodos/Desktop/Research/LINCS/data/expr/tensor/tsize/large/tensor_annot.RData')
L = GetLincsAnnot()
A = L[L$pert_id %in% annot$pertIds,]
targets_list = Str2Vec(na.omit(A$targets_entrez), split=':')
targets = unique(unlist(targets_list))

# tCount = list() #vector(mode='numeric', length=length(targets))
# for(target in targets){
#   tCount[[target]] = length(which(grepl(target, na.omit(A$targets_entrez))))
# }
names(targets) = MapEntrezToSymbol(targets, lm=FALSE)
tCount = sapply(targets, function(tar){length(which(sapply(targets_list, function(tlist) tar %in% tlist)))})
hist(tCount, breaks=50)

## Maybe should use SEA algorithm or the like to predict targets for the other
## guys (take all scores above certain stringent threshold)

# Then apply some network smoothing based on (gene-regulatory network?) to
# identify more connections between targets

atc = Str2Vec(na.omit(A$atc_code), split='[|]')
atc_letters = lapply(atc, function(ATC) sapply(ATC, function(y) substr(y, 1, 1)))

all_letters = sort(unique(unlist(atc_letters)))

n = list()
for(letter in all_letters){
  n[[letter]] = length(which(sapply(atc_letters, function(ATC) letter %in% ATC)))
}
