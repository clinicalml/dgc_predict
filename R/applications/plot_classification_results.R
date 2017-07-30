
# Plot results from run_binary_classifications.R

library(plotly)
library(ggplot2)
library(reshape2)
library(ggrepel)
library(plyr)

load('../results/classification/2017-07-30-03-47-00/results.RData')

Outcome2Category = function(outcome){
  return(sapply(as.character(outcome), function(x) unlist(strsplit(x, split='[.]'))[[1]]))
}

numSigs = c(MCF7=1505, VCAP=1368, PC3=1340, A375=1168, A549=1139, 
            HA1E=1127, HT29=1022, HCC515=934, HEPG2=798, NPC=441)


# Melt ROC results
R = melt(ROC)
names(R) = c('ROC', 'eval_type', 'model', 'obs_type', 'feature_type', 'outcome')
R$category = Outcome2Category(R$outcome)

# Subset to evaluations on measured signatures
Rmeas = RemoveDfColumns(subset(R, eval_type == 'eval_meas'), 'eval_type')
Rmeas$ROC[is.na(Rmeas$ROC)] = 0.5
Rmeas = split(Rmeas, Rmeas$obs_type)
Rmeas = merge(Rmeas$full, Rmeas$obs, by=c('feature_type', 'outcome','model'), 
              all=TRUE, suffixes=c('.full','.obs'))
Rmeas$delta = Rmeas$ROC.full - Rmeas$ROC.obs

for(threshold in c(0.5, 0.6, 0.7)){
  A = subset(Rmeas, ROC.full > threshold | ROC.obs > threshold)
  print(ggplot(A, aes(x=model, y=delta, group=model, fill=model)) + geom_boxplot() + 
          ggtitle(sprintf('Deltas per model and outcome, threshold = %0.1f', threshold)))
}

for(threshold in c(0.5, 0.6, 0.7)){
  A = subset(Rmeas, ROC.full > threshold | ROC.obs > threshold)
  print(ggplot(A, aes(x=reorder(feature_type, -delta, FUN=median), y=delta, group=feature_type, fill=feature_type)) + geom_boxplot() + 
          ggtitle(sprintf('Deltas per model and outcome, threshold = %0.1f', threshold))
        + theme(axis.text.x = element_text(angle = 45, hjust = 1)))
}

for(threshold in c(0.5, 0.6, 0.7)){
  A = subset(Rmeas, ROC.full > threshold | ROC.obs > threshold)
  print(ggplot(A, aes(x=outcome, y=delta, group=outcome, fill=category.obs)) + geom_boxplot() + 
          ggtitle(sprintf('Deltas per model and outcome, threshold = %0.1f', threshold))
        + theme(axis.text.x = element_text(angle = 45, hjust = 1)))
}

# Setup df so that I can color faceted plots by category
r = data.frame(outcome = unique(R$outcome))
r$category = Outcome2Category(r$outcome)
r$ROC.meas = r$ROC.imp = r$ROC.obs = r$ROC.full = 1 


# Plot AUCs on same evaluation set, either with or without using the predicted signatures
Rmeas = RemoveDfColumns(subset(R, eval_type == 'eval_meas'), 'eval_type')
Rmeas$ROC[is.na(Rmeas$ROC)] = 0.5
Rmeas = split(Rmeas, Rmeas$obs_type)
Rmeas = merge(Rmeas$full, Rmeas$obs, by=c('feature_type', 'outcome'), all=TRUE, suffixes=c('.full','.obs'))
p = ggplot(Rmeas, aes(x=ROC.obs, y=ROC.full)) +
  geom_rect(data=r, aes(fill=category), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf,
            alpha = 0.4) +
  geom_point() +
  geom_abline(intercept=0, slope=1) +
  geom_hline(yintercept=0.5, col='DarkGrey', linetype=2) +
  geom_vline(xintercept=0.5, col='DarkGrey', linetype=2) +
  geom_text_repel(data=Rmeas, aes(x=ROC.obs, y=ROC.full, label=feature_type), color='black',
                  size=2, segment.size = 0.3, segment.color='black', 
                  box.padding = unit(0.5, "lines"),
                  point.padding = unit(0.1, "lines")) +
  facet_wrap(~outcome) +
  scale_fill_brewer(palette='Set3') +
  xlab('AUC using ONLY MEASURED signatures') +
  ylab('AUC using MEASURED + IMPUTED signatures') +
  ggtitle('Does addition of imputed signatures improve classification accuracy?') +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlim(c(0,1)) + ylim(c(0,1))
print(p)
ggsave(plot=p, filename=paste0(resDir, '/Does_imputation_help.pdf'), width=15, height=12)

# # Make box plots of AUC deltas, per feature and per outcome
# # Filter by cases where neither made it above 0.5 (?)
# A = subset(R2, ROC.full > 0.6 | ROC.obs > 0.6)
# A$delta = A$ROC.full - A$ROC.obs
# ggplot(A, aes(x=reorder(outcome, delta, FUN=median), y=delta, group=outcome, fill=category.full)) +
#   geom_boxplot() +
#   ggtitle('Change in AUC per prediction task when adding imputed signatures')
# ggplot(A, aes(x=reorder(feature_type, delta, FUN=median), y=delta, group=feature_type))

# # Then see if predictions on imputed vs. measured signatures have comparable AUCs
# R3 = RemoveDfColumns(subset(R, obs_type == 'full'), 'obs_type')
# R3 = split(R3, R3$eval_type)
# R3 = merge(R3$eval_meas, R3$eval_imp, by=c('feature_type', 'outcome'), all=TRUE, suffixes=c('.meas','.imp'))
# p2 = ggplot(R3, aes(x=ROC.meas, y=ROC.imp)) + 
#   geom_rect(data=r, aes(fill=category), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf,
#             alpha = 0.4) +
#   geom_point() + 
#   geom_abline(intercept=0, slope=1) +
#   geom_hline(yintercept=0.5, col='DarkGrey', linetype=2) +
#   geom_vline(xintercept=0.5, col='DarkGrey', linetype=2) +
#   geom_text_repel(data=R3, aes(x=ROC.meas, y=ROC.imp, label=feature_type), color='black',
#                   size=2, segment.size = 0.3, segment.color='black', 
#                   box.padding = unit(0.5, "lines"),
#                   point.padding = unit(0.1, "lines")) +
#   facet_wrap(~outcome) +
#   scale_fill_brewer(palette='Set3') +
#   xlab('AUC evaluated on data points with MEASURED signatures') +
#   ylab('AUC evaluated on the remaining data points with IMPUTED signatures') +
#   ggtitle('Do predictions on imputed vs. measured signatures have comparable AUCs?') +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   xlim(c(0,1)) + ylim(c(0,1))
# print(p2)
# ggsave(plot=p2, filename=paste0(resDir, '/Is_accuracy_comparable.pdf'), width=15, height=12)
# 
# 
# ### Re-facet the first plot by feature type
# ggplot(R2, aes(x=ROC.obs, y=ROC.full)) +
#   #geom_rect(data=r, aes(fill=category), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf,
#   #          alpha = 0.4) +
#   geom_point() +
#   geom_abline(intercept=0, slope=1) +
#   geom_hline(yintercept=0.5, col='DarkGrey', linetype=2) +
#   geom_vline(xintercept=0.5, col='DarkGrey', linetype=2) +
#   geom_text_repel(data=R2, aes(x=ROC.obs, y=ROC.full, label=outcome), color='black',
#                   size=2, segment.size = 0.3, segment.color='black', 
#                   box.padding = unit(0.5, "lines"),
#                   point.padding = unit(0.1, "lines")) +
#   facet_wrap(~feature_type) +
#   scale_fill_brewer(palette='Set3') +
#   xlab('AUC using ONLY MEASURED signatures') +
#   ylab('AUC using MEASURED + IMPUTED signatures') +
#   ggtitle('Does addition of imputed signatures improve classification accuracy?') +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   xlim(c(0,1)) + ylim(c(0,1))
# #ggsave(filename=PlotDir())