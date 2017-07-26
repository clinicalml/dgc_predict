setwd('~/projects/dg/dgc_predict/')

### Setup parallelization and load dependencies
library(doParallel)
cl = makeCluster(26)
registerDoParallel(cl)

options(error=recover)
library(RSQLite)
source('R/src/DataProc.R')
source('R/src/Utils.R')

library(caret)
set.seed(123)

### Run params
debug = FALSE
save = FALSE
plot = FALSE
dim = 100

### Load gene expression features. For now, remove concatenated features since
### it takes longer, and remove cv features
print('loading data...')
load(DataDir('expr/drug/tensor_features_for_drug_property_prediction.RData'))
L = L[setdiff(names(L), 'cat')]
L = lapply(L, function(x) x[setdiff(names(x), 'cv')])


### Load label matrices
load(DataDir('all_labels/all_labels.RData'))

### Reduce both L and Y if in debug mode
if(debug){
  L = L[c(1,3)]
  for(i in 1:2){
    L[[i]]$obs = L[[i]]$obs[,1:5]
    L[[i]]$full = L[[i]]$full[,1:5]
  }
  
  Y = Y[,c(1, 19, 25)]
}

### Reduce dimensionality of input
L = lapply(L, function(xx) lapply(xx, function(x){


### Setup caret parameters
fitControl = trainControl(method=ifelse(debug, 'LOOCV', 'cv'), number=10, classProbs=TRUE, allowParallel=TRUE,
                          summaryFunction=twoClassSummary, savePredictions='final')
if(debug){
  pGrid = expand.grid(lambda=0, diagonal=TRUE)
}else{
  pGrid = expand.grid(lambda=seq(0, 1, 0.15), diagonal=c(TRUE, FALSE))
}


### For each input feature type, run SDA (shrinkage discriminant analysis)


#system.time(out <- foreach(y = colnames(Y), .packages=c('caret','ROCR')) %dopar% {
for(y in colnames(Y)){
  ROC = list()
  params = list()
  y_perts = rownames(Y)[!is.na(Y[,y])]

  for(f in names(L)){

    # Get measured drugs
    Xmeas = L[[f]]$obs[y_perts,]
    measPerts = rownames(Xmeas)[!is.na(Xmeas[,1])]
    impPerts = rownames(Xmeas)[is.na(Xmeas[,1])]
   
    for(subset in names(L[[f]])){
      print(sprintf('Predicting %s with %s (%s) features', y, f, subset))
      
      # Construct data matrix
      X = L[[f]][[subset]][y_perts,]
      X = X[!is.na(X[,1]),] # remove missing signatures in the case of only measured data

      # Identify side effect labels
      labels = factor(Y[rownames(X),y], levels=c(0,1), labels=c('neg','pos'))
      labels_v = as.numeric(labels) - 1
      names(labels_v) = rownames(X)
      
      # Run model
      capture.output(out <- caret::train(x=X, y=labels, preProcess='center', method='sda',
                                  trControl=fitControl, verbose=FALSE, tuneGrid=pGrid, metric='ROC'))
      rownames(out$pred) = rownames(X)[out$pred$rowIndex]
      
      roc = list(eval_all = max(out$results$ROC), eval_meas=ComputeAUC(est = out$pred[measPerts,'pos'], labels=labels_v[measPerts]))
      roc$eval_imp = ifelse(length(impPerts)>0, ComputeAUC(est=out$pred[impPerts,'pos'], labels=labels_v[impPerts]),NA)

      params[[y]][[f]][[subset]] = out$bestTune
      ROC[[y]][[f]][[subset]] = roc
      #params[[f]][[subset]] = out$bestTune
      #ROC[[f]][[subset]] = roc
    }
  }
  #return(list(p=params, R=ROC))
  #return(ROC)
}

if(debug){
  load(ResultsDir(sprintf('sda%s.RData', ifelse(debug, '_debug', ''))))
  stopifnot(identical(ROC_save, ROC))
  stopifnot(identical(params_save, params))
  print('Test passed')
}


if(save){
  save(ROC, params, file=ResultsDir(sprintf('sda%s.RData', ifelse(debug, '_debug', ''))))
}


### Let's plot these results
if(plot){
  R = melt(ROC)
  names(R) = c('ROC', 'eval_type', 'obs_type', 'feature_type', 'side_effect')
  R = split(R, R$obs_type)
  R = merge(R$full, R$obs, by=c('feature_type', 'side_effect'), all=TRUE, suffixes=c('.full','.obs'))
  ggplot(rr, aes(x=ROC.obs, y=ROC.full)) + geom_point() + geom_abline(intercept=0, slope=1)
}

