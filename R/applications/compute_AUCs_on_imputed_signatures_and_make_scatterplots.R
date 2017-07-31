# 
# load('../results/classification/2017-07-30-03-47-00/results_ROC_counts_params.RData')
# #save(ROC, counts, params, file='../results/classification/2017-07-30-03-47-00/results_ROC_counts_params.RData')
# R_old = melt(ROC)
# C_old = melt(counts)
# 
# ### Load gene expression features. For now, remove concatenated features since
# ### it takes longer, and remove cv features
# print('loading data...')
# load(DataDir('expr/drug/tensor_features_for_drug_property_prediction_10cells_knn.RData'))
# L = L[setdiff(names(L), c('allcell','pca200','pca978'))]
# 
# ### Load label matrices
# load(DataDir('all_labels/all_labels.RData'))
# 
# ROC = list()
# counts = list()
# subset = 'full'
# models = c('regLogistic', 'parRF', 'knn')
# 
# for(y in colnames(Y)){
#   y_perts = rownames(Y)[!is.na(Y[,y])]
#   
#   for(f in names(L)){
#     
#     # Among drugs with labels for output y, identify which ones have measured vs. only imputed signatures
#     Xmeas = na.omit(L[[f]]$obs)
#     meas = intersect(rownames(Xmeas), y_perts)
#     imp = setdiff(y_perts, rownames(Xmeas))
# 
#     # Construct data matrix
#     X = na.omit(L[[f]][[subset]][y_perts,])
#     
#     # And get corresponding labels (factor and numeric versions)
#     lab = factor(Y[rownames(X),y], levels=c(0,1), labels=c('neg','pos'))
#     lab_v = setNames(Y[rownames(X),y], rownames(X))
#     
#     # Get label counts
#     C = list(nPos_meas = length(which(lab_v[meas]==1)), nTot_meas = length(meas),
#              nPos_imp = length(which(lab_v[imp]==1)), nTot_imp = length(imp),
#              nPos = length(which(lab_v == 1)), nTot = length(lab_v))
#     counts[[y]][[f]][[subset]] = C
#     
#     # Run models
#     for(model in models){
#       if(length(unique(lab)) == 2 && C$nPos >= 5 && (C$nTot - C$nPos) >= 5){
#         print(sprintf('Computing AUC for %s with %s (%s) features and %s model', y, f, subset, model))
#         out = OUT[[y]][[f]][[subset]][[model]]
#         rownames(out$pred) = rownames(X)[out$pred$rowIndex]
#         roc = list(eval_meas = ComputeAUC(est = out$pred[meas,'pos'], lab=lab_v[meas], na.rm=TRUE),
#                    eval_imp=ComputeAUC(est=out$pred[imp,'pos'], lab=lab_v[imp], na.rm=TRUE))
#       }else{
#         roc = list(eval_meas=NA, eval_imp=NA)
#       }
#       ROC[[y]][[f]][[subset]][[model]] = roc
#     }
#   }
# }
# 
# 
# library(compare)
# 
# # Test ROC eval_meas
# R_new = subset(melt(ROC), L5 == 'eval_meas')
# R_old = subset(R_old, L3=='full' & L5 == 'eval_meas')
# stopifnot(compare(R_new, R_old, allowAll = TRUE)$result)
# 
# # Test counts eval_meas
# C_new = melt(counts)
# C_old = subset(C_old, L3 == 'full')
# stopifnot(compare(C_old, C_new, allowAll = TRUE)$result)
# 
# # Make scatter plots
# 
# Outcome2Category = function(outcome){
#   return(sapply(as.character(outcome), function(x) unlist(strsplit(x, split='[.]'))[[1]]))
# }
# 
# R = melt(ROC)
# R = dcast(R, L1 + L2 + L3 + L4 ~ L5)
# names(R) = c('outcome','feature','subset','model', 'AUC_imp', 'AUC_meas')
# R$category = Outcome2Category(R$outcome)
# R$outcome = sapply(R$outcome, function(x) unlist(strsplit(x, split='[.]'))[[2]])
# 
# R = ddply(R, c('outcome','feature','category','subset'), 
#           summarize, AUC_imp=median(AUC_imp), AUC_meas=median(AUC_meas))
# 
# # Merge R with counts
# C = C_new[,-3]
# names(C) = c('count','variable','feature','outcome')
# C = dcast(C, feature + outcome  ~ variable, value.var='count')
# C$outcome = sapply(C$outcome, function(x) unlist(strsplit(x, split='[.]'))[[2]])
# RC = merge(R, C, all=TRUE, by=c('outcome','feature'))
# #RC$positive_ratio = RC$nPos_imp / RC$nPos
# 
# RC = subset(RC, nPos_imp >= 3 & nPos_meas >= 3)
load(ResultsDir('classification/2017-07-30-03-47-00/RC.RData'))

RC$diff = RC$AUC_imp - RC$AUC_meas
ATC = subset(RC, category=='ATC' & outcome %in% c('L','C','D'))

ggplot(ATC, aes(x=AUC_meas, y=AUC_imp, color=outcome, label=feature)) + 
  geom_abline(slope=1, intercept=0, color='DarkGrey', lty='dashed') + geom_point(size=9, alpha=0.6) +
  xlim(c(0,1)) + ylim(c(0,1)) +
  geom_text(size=3, color='black', check_overlap=TRUE) + theme_bw() + 
  guides(color=guide_legend(title='ATC Code')) +
  xlab('AUC on measured signatures') + ylab('AUC on predicted signatures') + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=18),
        legend.text=element_text(size=12), legend.title=element_text(size=14))
ggsave(PlotDir('ATC_code_imputed_vs_measured.svg'), width=8, height=7)

Targets = subset(RC, category == 'Target')
ggplot(Targets, aes(x=AUC_meas, y=AUC_imp, color=outcome, label=feature)) + 
  geom_abline(slope=1, intercept=0, color='DarkGrey', lty='dashed') + geom_point(size=9, alpha=0.6) +
  xlim(c(0,1)) + ylim(c(0,1)) +
  geom_text(size=3, color='black', check_overlap=TRUE) + theme_bw() + guides(color=guide_legend(title='Target')) +
  xlab('AUC on measured signatures') + ylab('AUC on predicted signatures') + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=18),
        legend.text=element_text(size=12), legend.title=element_text(size=14))
ggsave(PlotDir('Target_imputed_vs_measured.svg'), width=8.5, height=7)


