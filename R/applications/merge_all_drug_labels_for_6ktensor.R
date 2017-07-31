
# targets
TG = read.table(DataDir('targets/top_targets_in_6ktensor.csv'), sep=',', header=TRUE, row.names=1)
colnames(TG) = paste0('Target.', colnames(TG))

# ATC codes
ATC = read.table(DataDir('atc/top_14_atc_codes_in_6k_tensor.csv'), sep=',', header=TRUE, row.names=1)
colnames(ATC) = paste0('ATC.', colnames(ATC))

# merge into a single matrix
Y = transform( merge(ATC, TG, by=0, all=TRUE), row.names=Row.names, Row.names=NULL)

# Save 
save(Y, file=DataDir('all_labels/ATC_and_targets_for_6k_tensor.RData'))
