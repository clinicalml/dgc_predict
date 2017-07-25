
testAll = FALSE
source('R/init.R')

library(sda)
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

### For each input feature type, run SDA (shrinkage discriminant analysis)
AUC = list()
acc = list()

for(s in 1:Ns){
  se = colnames(SE)[s]

  for(f in 1:5){#Nf){
    feat = names(X)[f]
    
    # Get measured drugs
    x = X[[f]]$obs[rownames(SE),]
    measuredDrugs = rownames(x)[!is.na(x[,1])]
    imputedDrugs = rownames(x)[is.na(x[,1])]
   
    for(t in 1:length(X[[f]])){
      t_type = names(X[[f]])[t]
      
      print(sprintf('Predicting %s with %s (%s) features', se, feat, t_type))
      
      # Construct data matrix
      x = X[[f]][[t]][rownames(SE),]
      x = x - apply(x, 2, mean, na.rm=TRUE) # center data
      x = x[!is.na(x[,1]),] # remove missing signatures in the case of only measured data

      # Identify side effect labels
      labels = as.factor(SE[rownames(x),s]+1)
      nPos = length(which(labels==1))
      nTot = length(labels)
      
      # Estimate optimal regularization
      reg = sda(x, labels)$regularization * 2 ### I'm doubling the estimated shrinkage parameters!!!

      # Perform leave-one-out CV
      pred = vector(length=nrow(x), mode='numeric')
      post = vector(length=nrow(x), mode='numeric')
      
      for(d in 1:nrow(x)){
        x_train = x[-d,]
        y_train = labels[-d]
        sda.fit = sda(x_train, y_train, lambda=reg['lambda'], lambda.var=reg['lambda.var'], 
                      lambda.freqs=reg['lambda.freqs'], verbose=FALSE)
        pred[d] = predict(sda.fit, x[d,,drop=FALSE], verbose=FALSE)$class 
        post[d] = predict(sda.fit, x[d,,drop=FALSE], verbose=FALSE)$posterior[,2]
      }
      
      names(pred) = rownames(x)
      names(post) = rownames(x)
      names(labels) = rownames(x)
      
      # Compute accuracy
      acc[[se]][[feat]][[t_type]]$measured = 1 - mean((pred[measuredDrugs] - as.numeric(labels[measuredDrugs]))^2)
      acc[[se]][[feat]][[t_type]]$imputed = 1 - mean((pred[imputedDrugs] - as.numeric(labels[imputedDrugs]))^2)
      acc[[se]][[feat]][[t_type]]$random = -666 ### NEED TO COMPUTE THIS
      
      # Compute AUC
      #AUC[[se]][[feat]][[t_type]]$measured = mean((pred[measuredDrugs] - as.numeric(labels[measuredDrugs]))^2)
      #AUC[[se]][[feat]][[t_type]]$imputed = mean((pred[imputedDrugs] - as.numeric(labels[imputedDrugs]))^2)
      
      print(sprintf('ACCURACY = %0.2f measured, %0.2f imputed (vs. random=%0.2f)', acc[[se]][[feat]][[t_type]]$measured,
                    acc[[se]][[feat]][[t_type]]$imputed))

    }
  }
}


