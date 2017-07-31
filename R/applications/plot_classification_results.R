
# Plot results from run_binary_classifications.R

library(plotly)
library(ggplot2)
library(reshape2)
library(ggrepel)
library(plyr)
library(ggpubr)

options(error=recover)
library(RSQLite)
source('R/src/DataProc.R')
source('R/src/Utils.R')

load('../results/classification/2017-07-30-03-47-00/results.RData')

numSigs = c(MCF7=1505, VCAP=1368, PC3=1340, A375=1168, A549=1139, 
            HA1E=1127, HT29=1022, HCC515=934, HEPG2=798, NPC=441)


# Melt AUC results
R = melt(ROC)
names(R) = c('AUC', 'eval', 'model', 'subset', 'feature', 'outcome')
R = SplitOutcome(R)

# Subset to top three represented ATCs
R = subset(R, outcome %in% c('L','C','D') | category == 'Target')

# Subset to evaluations on measured signatures and reformat data
Rmeas = RemoveDfColumns(subset(R, eval == 'eval_meas'), 'eval')
Rmeas$AUC[is.na(Rmeas$AUC)] = 0.5
Rmeas = dcast(Rmeas, model + feature + outcome + category ~ subset, value.var='AUC')
Rmeas = ChangeColumnName(Rmeas, from=c('full','obs'), to=c('AUC.full', 'AUC.obs'))
Rmeas$diff = Rmeas$AUC.full - Rmeas$AUC.obs

for(threshold in c(0.5, 0.6, 0.7)){
  A = subset(Rmeas, AUC.full > threshold | AUC.obs > threshold)
  print(ggplot(A, aes(x=model, y=diff, group=model, fill=model)) + geom_boxplot() + 
          ggtitle(sprintf('Deltas per model and outcome, threshold = %0.1f', threshold)) + stat_compare_means())
}

for(threshold in c(0.5, 0.6, 0.7)){
  A = subset(Rmeas, AUC.full > threshold | AUC.obs > threshold)
  print(ggplot(A, aes(x=reorder(feature, -diff, FUN=median), y=diff, group=feature, fill=feature)) + geom_boxplot() + 
          ggtitle(sprintf('Deltas per model and outcome, threshold = %0.1f', threshold))
        + theme(axis.text.x = element_text(angle = 45, hjust = 1)))
}

for(threshold in c(0.5, 0.6, 0.7)){
  A = subset(Rmeas, AUC.full > threshold | AUC.obs > threshold)
  print(ggplot(A, aes(x=outcome, y=diff, group=outcome, fill=category.obs)) + geom_boxplot() + 
          ggtitle(sprintf('Deltas per model and outcome, threshold = %0.1f', threshold))
        + theme(axis.text.x = element_text(angle = 45, hjust = 1)))
}

# Setup df so that I can color faceted plots by category
r = data.frame(outcome = unique(R$outcome))
r$category = Outcome2Category(r$outcome)
r$AUC.meas = r$AUC.imp = r$AUC.obs = r$AUC.full = 1 


# Plot AUCs on same evaluation set, either with or without using the predicted signatures
Rmeas = RemoveDfColumns(subset(R, eval == 'eval_meas'), 'eval')
Rmeas$AUC[is.na(Rmeas$AUC)] = 0.5
Rmeas = split(Rmeas, Rmeas$subset)
Rmeas = merge(Rmeas$full, Rmeas$obs, by=c('feature', 'outcome'), all=TRUE, suffixes=c('.full','.obs'))
p = ggplot(Rmeas, aes(x=AUC.obs, y=AUC.full)) +
  geom_rect(data=r, aes(fill=category), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf,
            alpha = 0.4) +
  geom_point() +
  geom_abline(intercept=0, slope=1) +
  geom_hline(yintercept=0.5, col='DarkGrey', linetype=2) +
  geom_vline(xintercept=0.5, col='DarkGrey', linetype=2) +
  geom_text_repel(data=Rmeas, aes(x=AUC.obs, y=AUC.full, label=feature), color='black',
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
# A = subset(R2, AUC.full > 0.6 | AUC.obs > 0.6)
# A$diff = A$AUC.full - A$AUC.obs
# ggplot(A, aes(x=reorder(outcome, diff, FUN=median), y=diff, group=outcome, fill=category.full)) +
#   geom_boxplot() +
#   ggtitle('Change in AUC per prediction task when adding imputed signatures')
# ggplot(A, aes(x=reorder(feature, diff, FUN=median), y=diff, group=feature))

# # Then see if predictions on imputed vs. measured signatures have comparable AUCs
# R3 = RemoveDfColumns(subset(R, subset == 'full'), 'subset')
# R3 = split(R3, R3$eval)
# R3 = merge(R3$eval_meas, R3$eval_imp, by=c('feature', 'outcome'), all=TRUE, suffixes=c('.meas','.imp'))
# p2 = ggplot(R3, aes(x=AUC.meas, y=AUC.imp)) + 
#   geom_rect(data=r, aes(fill=category), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf,
#             alpha = 0.4) +
#   geom_point() + 
#   geom_abline(intercept=0, slope=1) +
#   geom_hline(yintercept=0.5, col='DarkGrey', linetype=2) +
#   geom_vline(xintercept=0.5, col='DarkGrey', linetype=2) +
#   geom_text_repel(data=R3, aes(x=AUC.meas, y=AUC.imp, label=feature), color='black',
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
# ggplot(R2, aes(x=AUC.obs, y=AUC.full)) +
#   #geom_rect(data=r, aes(fill=category), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf,
#   #          alpha = 0.4) +
#   geom_point() +
#   geom_abline(intercept=0, slope=1) +
#   geom_hline(yintercept=0.5, col='DarkGrey', linetype=2) +
#   geom_vline(xintercept=0.5, col='DarkGrey', linetype=2) +
#   geom_text_repel(data=R2, aes(x=AUC.obs, y=AUC.full, label=outcome), color='black',
#                   size=2, segment.size = 0.3, segment.color='black', 
#                   box.padding = unit(0.5, "lines"),
#                   point.padding = unit(0.1, "lines")) +
#   facet_wrap(~feature) +
#   scale_fill_brewer(palette='Set3') +
#   xlab('AUC using ONLY MEASURED signatures') +
#   ylab('AUC using MEASURED + IMPUTED signatures') +
#   ggtitle('Does addition of imputed signatures improve classification accuracy?') +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   xlim(c(0,1)) + ylim(c(0,1))
# #ggsave(filename=PlotDir())