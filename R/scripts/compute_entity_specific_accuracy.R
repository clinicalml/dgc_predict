library(ggplot2)
library(dplyr)
library(plotrix)
library(R.matlab)
library(tidyr)

#### Parameters for DEG computation
normGene = FALSE
degMethod = 'sig'
perc = 2
symmetric = FALSE

#### Run parameters
saveToFile = TRUE
dim_compute = c('gene', 'drug', 'cell')
computeAUC = FALSE
computeCor = TRUE

#### LOAD DATA 
n = lapply(dimnames(tensors$meas), function(x) length(x))
imputation_methods = c('mean', 'mean2', 'knn', 'tensor')

##### COMPUTE AUC'S AND/OR CORRELATIONS FOR EACH ELEMENT OF EACH DIMENSION 

if(computeAUC){
  auc = list()
  for(method in imputation_methods){
    auc[[method]] = list()
  }
}

if(computeCor){
  C = list()
  for(method in imputation_methods){
    C[[method]] = list()
  }
}

numSig = list()
numDEG = list()
percDEG = list()

X = tensors$cv # list of cross-validated tensors
x = list() # vectorized versions of tensors

# compute labels
X$DEG = abs(TensorDEG(tensors$meas, normGene=normGene, method=degMethod, percDEG=perc, symmetric=symmetric)$D)
X$true = tensors$meas

for(dim in dim_compute){

  # permute matrices so that chosen dimension is first
  for(tensor in names(X)){
    X[[tensor]] = aperm(X[[tensor]], c(dim, setdiff(c('drug', 'gene','cell'), dim)))
  }

  for(idx in 1:n[[dim]]){
    print(idx)
    for(tensor in names(X)){
      x[[tensor]] = as.vector(X[[tensor]][idx,,])
    }

    num_entries = length(which(!is.na(x$true)))
    if(dim %in% c('drug', 'cell')){
      num_sig = num_entries / n$gene
    }else{
      num_sig = num_entries
    }
    numSig[[dim]][idx] = num_sig
    numDEG[[dim]][idx] = length(which(x$DEG==1))
    percDEG[[dim]][idx] = length(which(x$DEG==1)) / length(which(!is.na(x$DEG)))
    
    for(method in imputation_methods){
      if(computeAUC){
        if(all(c(0,1) %in% unique(x$DEG))){
          auc[[method]][[dim]][idx] = ComputeAUC(x[[method]], x$DEG)
        }else{
          auc[[method]][[dim]][idx] = NA
        }
      }
      
      if(computeCor){
        C[[method]][[dim]][idx] = cor(x[[method]], x$true, use='pairwise.complete')
      }
    }
  }
}

# and save to file
if(saveToFile){
   save(C, file=DataDir('results/tsize/large/accuracy_per_mode.RData'))
}
