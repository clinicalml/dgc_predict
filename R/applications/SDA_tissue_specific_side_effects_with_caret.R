
testAll = FALSE
source('R/init.R')

#library(sda)
library(caret)
set.seed(123)

### Load gene expression features
print('loading data...')
load(DataDir('expr/drug/tensor_features_for_drug_property_prediction.RData'))

### For now, remove concatenated features since it takes longer
X = X[setdiff(names(X), 'cat')]
Nf = length(X)

### Also remove cv features
X = lapply(X, function(x) x[setdiff(names(x), 'cv')])

### Load SE matrices
SE = read.table(DataDir('side_effect/sider_sliced.csv'), sep=',', header=TRUE, row.names=1)
Ns = ncol(SE)

### TODO: Concatenate side effects, ADRs, targets, and ATC codes together

### Setup caret parameters
fitControl = trainControl(method = 'cv', number = 10, classProbs = TRUE,
                          summaryFunction = twoClassSummary, savePredictions='final')
pGrid = expand.grid(lambda=seq(0, 1, 0.15), diagonal=c(TRUE, FALSE))

### For each input feature type, run SDA (shrinkage discriminant analysis)
ROC = list()
params = list()
#nG = 5

for(s in 1:Ns){
  se = colnames(SE)[s]
  ROC[[se]] = list()
  
  for(f in 1:Nf){
    
    feat = names(X)[f]
    ROC[[se]][[feat]] = list()
    
    # Get measured drugs
    x = X[[f]]$obs[rownames(SE),]
    measIdx = which(!is.na(x[,1]))
    impIdx = which(is.na(x[,1]))
   
    for(t in 1:length(X[[f]])){
      t_type = names(X[[f]])[t]
      ROC[[se]][[feat]][[t_type]] = list()
      
      print(sprintf('Predicting %s with %s (%s) features', se, feat, t_type))
      
      # Construct data matrix
      x = X[[f]][[t]][rownames(SE),]
      x = x[!is.na(x[,1]),] # remove missing signatures in the case of only measured data

      # Identify side effect labels
      labels = factor(SE[rownames(x),s], levels=c(0,1), labels=c('neg','pos'))
      
      # Run model
      capture.output(out <- train(x=x, y=labels, preProcess='center', method='sda', trControl=fitControl,
                  verbose=FALSE, tuneGrid=pGrid, metric='Kappa'))
      
      cvMeas = which(out$pred$rowIndex %in% measIdx)
      cvImp = which(out$pred$rowIndex %in% impIdx)
      
      ROC[[se]][[feat]][[t_type]]$all = max(out$results$ROC)
      ROC[[se]][[feat]][[t_type]]$meas = ComputeAUC(est = out$pred[cvMeas,'pos'], 
                                                    labels = as.numeric(out$pred$obs[cvMeas]) - 1)
      
      if(length(cvImp) > 0){
        ROC[[se]][[feat]][[t_type]]$imp = ComputeAUC(est = out$pred[cvImp,'pos'],
                                                     labels = as.numeric(out$pred$obs[cvImp]) - 1)
      }else{
        ROC[[se]][[feat]][[t_type]]$imp = NA
      }

      params[[se]][[feat]][[t_type]] = out$bestTune
      
    }
  }
}

save(ROC, params, file=ResultsDir('applications/sda_side_effects_v2.RData'))

### Let's plot these results
R = melt(ROC)
names(R) = c('ROC', 'eval', 'obs', 'feature', 'label')


Rmeas = subset(R, eval=='eval_meas')
Rmeas$obs

S = cast(R, ROC ~ obs)
#R = split(R, R$obs)
#R = merge(R$full, R$obs, by=c('feature_type', 'side_effect'), all=TRUE, suffixes=c('.full','.obs'))
p = ggplot(R, aes(x=ROC.obs, y=ROC.full)) + geom_point() + geom_abline(intercept=0, slope=1)
p = p + facet_wrap()
p <- p + facet_wrap( ~ day, ncol=2)
