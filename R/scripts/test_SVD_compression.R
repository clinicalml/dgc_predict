
dir = '/Users/rhodos/Desktop/Research/TC2_Lincs/data/ManyDrug/P0.1/matrix/onlyDrugs/'
fileTrain = paste0(dir, 'dataTrain_sm.tsv')
fileTest = paste0(dir, 'dataValidation_sm.tsv')
X = as.matrix(read.table(fileTrain, header=FALSE, sep='\t'))
Xtest = as.matrix(read.table(fileTest, header=FALSE, sep='\t'))

maxPC = 20
out = svd(X, nu=maxPC, nv=maxPC)

d10 = out$d
d10[11:length(d10)] = 0
D = diag(d10)
X10 = out$u %*% D[1:maxPC,1:maxPC] %*% t(out$v)
X10_test = Xtest %*% out$v[,1:10] %*% t(out$v[,1:10])

d20 = out$d
d20[21:length(d20)] = 0
D = diag(d20)
X20 = out$u %*% D[1:maxPC, 1:maxPC] %*% t(out$v)
X20_test = Xtest %*% out$v %*% t(out$v)

print(sprintf('Training error at rank 10 is %0.3f', ComputeSqrtErrRate(X, X10)))
print(sprintf('Training error at rank 20 is %0.3f', ComputeSqrtErrRate(X, X20)))
print(sprintf('Test error at rank 10 is %0.3f', ComputeSqrtErrRate(Xtest, X10_test)))
print(sprintf('Test error at rank 20 is %0.3f', ComputeSqrtErrRate(Xtest, X20_test)))
