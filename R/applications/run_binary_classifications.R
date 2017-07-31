setwd('~/projects/dg/dgc_predict/')

# ### Setup parallelization and load dependencies
library(doParallel)
cl = makeCluster(20)
registerDoParallel(cl)

options(error=recover)
library(RSQLite)
source('R/src/DataProc.R')
source('R/src/Utils.R')

library(caret)
set.seed(123)

### Run params
debug = FALSE
save = TRUE
models = c('regLogistic', 'parRF', 'knn', 'sda')

### Setup results dir
if(!debug){
  resDir = ResultsDir(paste0('classification/', gsub(Sys.time(), pattern=':| ', replacement='-')))
  MakeDir(resDir)
  from=paste0(getwd(), '/R/applications/run_binary_classifications.R')
  stopifnot(file.copy(from=from, to=resDir, overwrite=TRUE))
}else{
  resDir = ResultsDir(paste0('classification/debug'))
}

### Load gene expression features. For now, remove concatenated features since
### it takes longer, and remove cv features
print('loading data...')
load(DataDir('expr/drug/tensor_features_for_drug_property_prediction_10cells_knn.RData'))
L = L[setdiff(names(L), c('allcell','pca200','pca978'))]

### Load label matrices
load(DataDir('all_labels/all_labels.RData'))

### Reduce both L and Y if in debug mode
if(debug){
  L = L[c('mean','MCF7')]
  for(i in 1:2){
     L[[i]]$obs = L[[i]]$obs[,1:5]
     L[[i]]$full = L[[i]]$full[,1:5]
  }
  #Y = Y[,c('Target.CYP2C19','ATC.L')]
  Y = Y[,c(1, 19, 25)]
}

### Setup caret parameters
fitControl = trainControl(number=10, repeats=1, classProbs=TRUE, allowParallel=TRUE,
                          summaryFunction=twoClassSummary, savePredictions='final', 
                          method='cv')#ifelse(debug, 'LOOCV', 'cv'))

pGrid = list()
if(debug){
  pGrid$sda = expand.grid(lambda=0, diagonal=TRUE)
  pGrid$regLogistic = expand.grid(cost = 2^15, loss = "L2_dual", epsilon = 0.1)
  pGrid$parRF = expand.grid(mtry=4)
  pGrid$knn = expand.grid(k=1)
}else{
  pGrid$sda = expand.grid(lambda=seq(0, 1, 0.15), diagonal=c(TRUE, FALSE))
  pGrid$regLogistic = expand.grid(cost = 2^(-10:20), loss = "L2_dual", epsilon = 0.1)
  pGrid$parRF = expand.grid(mtry=seq(10, 40, 10))
  pGrid$knn = expand.grid(k=seq(1:10))
}

###########################  RUN ################################
ROC = list()
params = list()
counts = list()
OUT = list()

for(y in colnames(Y)){
  y_perts = rownames(Y)[!is.na(Y[,y])]

  for(f in names(L)){

    # Among drugs with labels for output y, identify which ones have measured vs. only imputed signatures
    Xmeas = na.omit(L[[f]]$obs)
    meas = intersect(rownames(Xmeas), y_perts)
    imp = setdiff(y_perts, rownames(Xmeas))
   
    for(subset in c('obs','full')){
      
      # Construct data matrix
      X = na.omit(L[[f]][[subset]][y_perts,])

      # And get corresponding labels (factor and numeric versions)
      lab = factor(Y[rownames(X),y], levels=c(0,1), labels=c('neg','pos'))
      lab_v = setNames(Y[rownames(X),y], rownames(X))
      
      # Get label counts
      C = list(nPos_meas = length(which(lab_v[meas]==1)), nTot_meas = length(meas),
               nPos_imp = length(which(lab_v[imp]==1)), nTot_imp = length(imp),
               nPos = length(which(lab_v == 1)), nTot = length(lab_v))
      counts[[y]][[f]][[subset]] = C

      # Run models
      for(model in models){
        if(length(unique(lab)) == 2 && C$nPos >= 5 && (C$nTot - C$nPos) >= 5){
          print(sprintf('Predicting %s with %s (%s) features and %s model', y, f, subset, model))
          capture.output(out <- caret::train(x=X, y=lab, preProcess='center', method=model, trControl=fitControl, tuneGrid=pGrid[[model]], metric='ROC'))
          rownames(out$pred) = rownames(X)[out$pred$rowIndex]
          
           roc = list(eval_all = ComputeAUC(est = out$pred[names(lab_v),'pos'], lab=lab_v, na.rm=TRUE),
                     eval_meas = ComputeAUC(est = out$pred[meas,'pos'], lab=lab_v[meas], na.rm=TRUE))
           
           # if(subset == 'full'){
           #   X_cv = L[[f]]$cv[meas,]
           #   pred = predict(out, X_cv, type = "prob")
           #   roc$eval_imp = ComputeAUC(est=pred$pos, lab=lab_v[meas], na.rm=TRUE)
           #   print(sprintf('  ROC = %0.2f (imp)', roc$eval_imp))
           # }
       
           prms = out$bestTune
           OUT[[y]][[f]][[subset]][[model]] = out
           print(sprintf('  ROC = %0.2f (all), %0.2f (meas)', roc$eval_all, roc$eval_meas))
        }else{
          roc = list(eval_all=NA, eval_meas=NA, eval_imp=NA)
          prms = suppressWarnings(pGrid[[model]][1,] + NA)
        }
        params[[y]][[f]][[subset]][[model]] = prms
        ROC[[y]][[f]][[subset]][[model]] = roc
      }
    }
  }
  
  print('saving...')
  if(save){save(ROC, params, counts, OUT, file=paste0(resDir, '/results.RData'))}
  print('...done!')
}

###########################  TEST  ################################

# if(debug){
#   load(sprintf('%s/sda%s.RData', resDir, ifelse(debug, '_debug', '')))
#   stopifnot(unlist(compare::compare(ROC, ROC_save, allowAll = TRUE)$result))
#   stopifnot(unlist(compare::compare(params, params_save, allowAll = TRUE)$result))
#   print('Test passed')
# }


### Copy this script to the results dir
from=paste0(getwd(), '/R/applications/run_binary_classifications.R')
stopifnot(file.copy(from=from, to=resDir, overwrite=TRUE))
