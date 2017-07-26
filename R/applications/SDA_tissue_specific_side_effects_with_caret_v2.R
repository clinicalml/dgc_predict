setwd('~/projects/dg/dgc_predict/')

### Setup parallelization and load dependencies
library(doParallel)
cl = makeCluster(10)
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
#L = lapply(L, function(xx) lapply(xx, function(x){


### Setup caret parameters
fitControl = trainControl(method=ifelse(debug, 'LOOCV', 'cv'), number=10, classProbs=TRUE, allowParallel=TRUE,
                          summaryFunction=twoClassSummary, savePredictions='final')
if(debug){
  pGrid = expand.grid(lambda=0, diagonal=TRUE)
}else{
  pGrid = expand.grid(lambda=seq(0, 1, 0.15), diagonal=c(TRUE, FALSE))
}


### For each input feature type, run model
ROC = list()
params = list()

for(y in colnames(Y)){
  y_perts = rownames(Y)[!is.na(Y[,y])]

  for(f in names(L)){

    # Among drugs with labels for output y, identify which ones have measured vs. only imputed signatures
    Xmeas = L[[f]]$obs[y_perts,]
    meas = rownames(Xmeas)[!is.na(Xmeas[,1])]
    imp = rownames(Xmeas)[is.na(Xmeas[,1])]
   
    for(subset in names(L[[f]])){
      print(sprintf('Predicting %s with %s (%s) features', y, f, subset))
      
      # Construct data matrix
      X = na.omit(L[[f]][[subset]][y_perts,])

      # And get corresponding labels (factor and numeric versions)
      lab = factor(Y[rownames(X),y], levels=c(0,1), labels=c('neg','pos'))
      lab_v = setNames(Y[rownames(X),y], rownames(X))

      # Run model
      if(length(unique(lab)) == 2){
        capture.output(out <- caret::train(x=X, y=lab, preProcess='center', method='sda', trControl=fitControl, tuneGrid=pGrid, metric='ROC'))
        rownames(out$pred) = rownames(X)[out$pred$rowIndex]
        roc = list(eval_all = ComputeAUC(est = out$pred[names(lab_v),'pos'], lab=lab_v, na.rm=TRUE),
                   eval_meas = ComputeAUC(est = out$pred[meas,'pos'], lab=lab_v[meas], na.rm=TRUE),
                   eval_imp = ifelse(length(na.omit(out$pred[imp,'pos'])>0), ComputeAUC(est=out$pred[imp,'pos'], lab=lab_v[imp], na.rm=TRUE),NA))
      }else{
        roc = list(eval_all=NA, eval_meas=NA, eval_imp=NA)
      }

      params[[y]][[f]][[subset]] = out$bestTune
      ROC[[y]][[f]][[subset]] = roc
    }
  }
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
Outcome2Category = function(outcome){
  out = sapply(outcome, function(x) unlist(strsplit(x, split='[.]'))[[1]])
  return(out)
}


### Setup for plotting
library(plotly)
library(ggplot2)
library(reshape2)

R = melt(ROC)
names(R) = c('ROC', 'eval_type', 'obs_type', 'feature_type', 'outcome')
R$category = Outcome2Category(R$outcome)

r = data.frame(outcome = unique(R$outcome))
r$category = Outcome2Category(r$outcome)
r$ROC.meas = r$ROC.imp = r$ROC.obs = r$ROC.full = 1 

# Plot AUCs on same evaluation set, either with or without using the predicted signatures
R2 = RemoveDfColumns(subset(R, eval_type == 'eval_meas'), 'eval_type')
R2 = split(R2, R2$obs_type)
R2 = merge(R2$full, R2$obs, by=c('feature_type', 'outcome'), all=TRUE, suffixes=c('.full','.obs'))
p = ggplot(R2, aes(x=ROC.obs, y=ROC.full)) +
  geom_rect(data=r, aes(fill=category), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf,
          alpha = 0.4) +
  geom_point() + 
  geom_abline(intercept=0, slope=1) + 
  geom_hline(yintercept=0.5, col='DarkGrey', linetype=2) +
  geom_vline(xintercept=0.5, col='DarkGrey', linetype=2) +
  facet_wrap(~outcome) + 
  scale_fill_brewer(palette='Set3') +
  xlab('AUC using ONLY MEASURED signatures') +
  ylab('AUC using MEASURED + IMPUTED signatures') +
  ggtitle('Does addition of imputed signatures improve classification accuracy?') +
  theme(plot.title = element_text(hjust = 0.5))
ggplotly(p)


# Then see if predictions on imputed vs. measured signatures have comparable AUCs
R3 = RemoveDfColumns(subset(R, obs_type == 'full'), 'obs_type')
R3 = split(R3, R3$eval_type)
R3 = merge(R3$eval_meas, R3$eval_imp, by=c('feature_type', 'outcome'), all=TRUE, suffixes=c('.meas','.imp'))
p2 = ggplot(R3, aes(x=ROC.meas, y=ROC.imp)) + 
  geom_rect(data=r, aes(fill=category), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf,
            alpha = 0.4) +
  geom_point() + 
  geom_abline(intercept=0, slope=1) +
  geom_hline(yintercept=0.5, col='DarkGrey', linetype=2) +
  geom_vline(xintercept=0.5, col='DarkGrey', linetype=2) +
  facet_wrap(~outcome) +
  scale_fill_brewer(palette='Set3') +
  xlab('AUC evaluated on data points with MEASURED signatures') +
  ylab('AUC evaluated on the remaining data points with IMPUTED signatures') +
  ggtitle('Do predictions on imputed vs. measured signatures have comparable AUCs?') +
  theme(plot.title = element_text(hjust = 0.5))
ggplotly(p2)
