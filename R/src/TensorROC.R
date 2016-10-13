library(ggplot2)
#library(plyr)

# choose indices between 1 and n_long so that I have approximately n_short
# indices evenly spaced and including the first and last points
SubsetIdx <- function(n_long, n_short){
  n_skip = floor(n_long / n_short)
  idx = seq(1, n_long, n_skip)
  if( idx[length(idx)] != n_long){
    idx = c(idx, n_long)
  }
  return(idx)
}

SubsetRoc <- function(roc.in, n_short){
  if(length(roc.in) > 2*n_short){
    n_long = length(roc.in$tpr)
    idx = SubsetIdx(n_long, n_short)
    roc.out = list()
    roc.out$tpr = roc.in$tpr[idx]
    roc.out$fpr = roc.in$fpr[idx]
  }else{
    roc.out = roc.in
  }
  return(roc.out)
}

# RunGseaV2 <- function(profile, k=50, up=TRUE){
#   if(up){
#     genes = names(sort(profile, decreasing=TRUE)[1:k])
#   }else{
#     genes = names(sort(profile)[1:k])
#   }
#   gsea = RunGsea(profile=profile, dz_genes_up=genes, dz_genes_down=c(),
#                  gs='go_kegg', minGeneSetSize=2, doGSOA=FALSE, doGSEA=TRUE)
#   #out = GetGSEAResultsV2(gsea)
#   return(gsea)
# }

