load('../results/classification/2017-07-30-03-47-00/new_ROC_and_counts_with_correct_AUC_imp_values.RData')

Outcome2Category = function(outcome){
  return(sapply(as.character(outcome), function(x) unlist(strsplit(x, split='[.]'))[[1]]))
}

# Reformat results
R = melt(ROC)
R = dcast(R, L1 + L2 + L3 + L4 ~ L5)
names(R) = c('outcome','feature','subset','model', 'AUC_imp', 'AUC_meas')
R$category = Outcome2Category(R$outcome)
R$outcome = sapply(R$outcome, function(x) unlist(strsplit(x, split='[.]'))[[2]])

# Merge R with counts
C = melt(counts)[,-3]
names(C) = c('count','variable','feature','outcome')
C = dcast(C, feature + outcome  ~ variable, value.var='count')
C$outcome = sapply(C$outcome, function(x) unlist(strsplit(x, split='[.]'))[[2]])
RC = merge(R, C, all=TRUE, by=c('outcome','feature'))
RC = subset(RC, category=='Target' | outcome %in% c('L','C','D')) 
print(sprintf('Starting with %d experiments', nrow(RC)))

# Filter by two different thresholds to avoid 'non-signal' cases
RC = subset(RC, AUC_imp > 0.5 | AUC_meas > 0.5 )
print(sprintf('After filtering by AUC > 0.5, %d experiments remaining', nrow(RC)))

RC = subset(RC, nPos_imp >= 3 & nPos_meas >= 3)
print(sprintf('After filtering by num labels >= 3, %d experiments remaining', nrow(RC)))

Targets = subset(RC, category == 'Target')
ATC = subset(RC, category == 'ATC')
save(RC, ATC, Targets, file=ResultsDir('classification/2017-07-30-03-47-00/RC_all_models.RData'))

t.test(ATC$AUC_imp, ATC$AUC_meas, paired=TRUE)
t.test(Targets$AUC_imp, Targets$AUC_meas, paired=TRUE)
