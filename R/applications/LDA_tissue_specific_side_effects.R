
testAll = FALSE
source('R/init.R')

library(GeoDE)
set.seed(123)

### Load gene expression features
print('loading data...')
load(DataDir('expr/drug/tensor_features_for_drug_property_prediction.RData'))

### For now, remove concatenated features since it takes longer
X = X[setdiff(names(X), 'cat')]

### Load SE matrices
#AE = read.table(DataDir('side_effect/faers_sliced.csv'), sep=',', header=TRUE, row.names=1)
SE = read.table(DataDir('side_effect/sider_sliced.csv'), sep=',', header=TRUE, row.names=1)

### TODO: Concatenate side effects, ADRs, targets, and ATC codes together

### For each input feature type, run characteristic direction 
out = list()
for(s in 1:ncol(SE)){
  se = colnames(SE)[s]

  for(f in 1:length(X)){
    feat = names(X)[f]
   
    for(t in 1:length(X[[f]])){
      t_type = names(X[[f]])[t]
      
      print(sprintf('Predicting %s with %s (%s) features', se, feat, t_type))
      
      # Construct data matrix
      x = X[[f]][[t]][rownames(SE),]
      x = x - apply(x, 2, mean, na.rm=TRUE) # center data
      x = x[!is.na(x[,1]),] # remove missing signatures in the case of only measured data
      data = data.frame(genenames=colnames(x), t(x))
      colnames(data) = gsub(colnames(data), pattern='[.]', replacement='-')
      
      # Identify side effect labels
      labels = as.factor(SE[colnames(data[,-1]),s] + 1)
      nPos = length(which(labels==2))
      nTot = length(labels)
      
      # Run characteristic direction
      w = chdirAnalysis(data, labels, CalculateSig=FALSE)$chdirprops$chdir[[1]]
      
      # Project data onto characteristic direction vector
      proj_x = t(w) %*% t(x)
      
      # Plot
      y = ifelse(labels==1, 1, 1.05)
      plot(x=proj_x, y=y, col=ifelse(labels==1, 'blue', 'red'), pch=8, ylim=c(0,2),
           main=sprintf('%s (%d/%d), using %s (%s) input', se, nPos, nTot, feat, t_type))
      
      # Run LDA
      #out[[se]][[feat]] = lda(x=x,grouping=labels, CV=TRUE)
      # Error: variables are collinear
      
    }
  }
}


