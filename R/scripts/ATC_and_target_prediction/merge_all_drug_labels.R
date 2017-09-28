load(DataDir('expr/tensor_features_for_drug_property_prediction_10cells_knn.RData'))

# targets
TG = read.table(DataDir('targets/top_7_targets_in_large_tensor.csv'), sep=',', header=TRUE, row.names=1)
colnames(TG) = paste0('Target.', colnames(TG))

# ATC codes
ATC = read.table(DataDir('atc/top_6_atc_codes_in_large_tensor.csv'), sep=',', header=TRUE, row.names=1)
colnames(ATC) = paste0('ATC.', colnames(ATC))

# merge into a single matrix
Y = transform( merge(ATC, TG, by=0, all=TRUE), row.names=Row.names, Row.names=NULL)
Y = Y[intersect(rownames(Y), rownames(L[[1]]$full)),]

save(Y, file=DataDir('all_labels/all_labels.RData'))
