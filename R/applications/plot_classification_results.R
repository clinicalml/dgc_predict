
# Plot results from run_binary_classifications.R

library(plotly)
library(ggplot2)
library(reshape2)
library(ggrepel)
library(plyr)
library(ggpubr)
library(RColorBrewer)

options(error=recover)
library(RSQLite)
source('R/src/DataProc.R')
source('R/src/Utils.R')

load('../results/classification/2017-07-30-03-47-00/results_ROC_counts_params.RData')

SplitOutcome = function(df){
  df$category = sapply(as.character(df$outcome), function(x) unlist(strsplit(x, split='[.]'))[[1]])
  df$outcome = sapply(as.character(df$outcome), function(x) unlist(strsplit(x, split='[.]'))[[2]])
  return(df)
}

FilterByAUC = function(df, threshold=0.5){
  return(subset(df, AUC.full > threshold | AUC.obs > threshold))
}

threshold = 0.5 # results where AUCs were both below this in the two comparisons are thrown out

# Melt AUC results
R = melt(ROC)
names(R) = c('AUC', 'eval', 'model', 'subset', 'feature', 'outcome')

# Setup df so that I can color faceted plots by category
r = SplitOutcome(data.frame(outcome = unique(R$outcome)))
r$AUC.meas = r$AUC.imp = r$AUC.obs = r$AUC.full = 1
r = subset(r, outcome %in% c('L','C','D') | category == 'Target')

# Split outcome into outcome and category
R = SplitOutcome(R)

# Subset to top three represented ATCs
R = subset(R, outcome %in% c('L','C','D') | category == 'Target')

# Subset to evaluations on measured signatures and reformat data
Rmeas = RemoveDfColumns(subset(R, eval == 'eval_meas'), 'eval')
Rmeas$AUC[is.na(Rmeas$AUC)] = 0.5
Rmeas = dcast(Rmeas, model + feature + outcome + category ~ subset, value.var='AUC')
Rmeas = ChangeColumnName(Rmeas, from=c('full','obs'), to=c('AUC.full', 'AUC.obs'))
Rmeas$diff = Rmeas$AUC.full - Rmeas$AUC.obs

# Add numSigs
numSigs = c(MCF7=1505, VCAP=1368, PC3=1340, A375=1168, A549=1139, 
            HA1E=1127, HT29=1022, HCC515=934, HEPG2=798, NPC=441)
Rmeas$num_sigs = numSigs[Rmeas$feature]

### First compute overall difference between the two groups
subset(Rmeas, AUC.full > threshold )

### Plot results per model
p = list()
for(threshold in thresholds){
  A = subset(Rmeas, AUC.full > threshold | AUC.obs > threshold)
  print(sprintf('%d of %d kept at threshold=%0.2f', nrow(A), nrow(Rmeas), threshold))
  print(ggplot(A, aes(x=model, y=diff, group=model, fill=model)) + geom_boxplot() +
          ggtitle(sprintf('Deltas per model, threshold = %0.1f', threshold)) + stat_compare_means())
  ggsave(PlotDir(sprintf('Deltas_per_model_threshold_%0.2f.pdf', threshold)))
  p[[as.character(threshold)]] = lapply(split(A, A$model), function(x) t.test(x$diff)$p.value)
}

### Plot results per feature
p_feature = p_outcome = list()
idx = Rmeas$category == 'ATC'
Rmeas$outcome[idx] = paste0('ATC ' , Rmeas$outcome[idx])
for(threshold in thresholds[1]){
  A = subset(Rmeas, AUC.full > threshold | AUC.obs > threshold)
  print(sprintf('%d of %d kept at threshold=%0.2f', nrow(A), nrow(Rmeas), threshold))
  p1 = ggplot(A, aes(x=reorder(feature, -diff, FUN=median), y=diff, group=feature, fill=num_sigs)) +
          geom_hline(yintercept = 0, color='grey', lwd=1) +
          geom_boxplot(alpha=0.9) +
          ggtitle(sprintf('Deltas per feature, threshold = %0.1f', threshold)) +
          theme_bw() + ylim(c(-0.5,0.5)) + scale_fill_gradientn(colours=brewer.pal(5,'YlOrRd')) +
          theme(text = element_text(size=18), axis.text.x = element_text(angle = 45, hjust = 1),
                plot.title = element_text(size=20, hjust=0.5),
                legend.justification='bottom', legend.position='bottom') + 
               # legend.justification=c(0.05,0.05), legend.position=c(0.05,0.05)) + # legend.position='bottom') +
          labs(fill='# Signatures  \n  Measured\n', x='', y=expression(AUC [impute] - AUC [orig]),
                title='Improvement in AUC per Drug Signature Type')
  #ggsave(PlotDir('DeltaAUC_per_Feature.svg'), height=7, width=8)
  # p_feature[[as.character(threshold)]] = lapply(split(A, A$feature), function(x){
  #   if(nrow(x) > 1){out = t.test(x$diff)$p.value}
  #   else{out = NA}
  #   return(out)
  # })
  
  p2 =  ggplot(A, aes(x=reorder(outcome, -diff, FUN=median), y=diff, group=outcome, fill=category)) + ylim(c(-0.5,0.5)) +
      geom_hline(yintercept = 0, color='grey', lwd=1) + geom_boxplot(alpha=0.8) + scale_fill_manual(values=c('#CC9966','#009999')) +
           ggtitle(sprintf('Deltas per outcome, threshold = %0.1f', threshold)) + theme_bw() +
      labs(x='', y=expression(AUC [impute] - AUC [orig]), title='Improvement in AUC per Prediction Task',
           fill='Category') +
      theme(text = element_text(size=18), axis.text.x = element_text(angle = 45, hjust = 1), 
            plot.title = element_text(size=20, hjust=0.5),
            legend.justification='bottom', legend.position='bottom')
  #ggsave(PlotDir('DeltaAUC_per_Outcome.svg'), height=7, width=8)
  # p_outcome[[as.character(threshold)]] = lapply(split(A, A$outcome), function(x){
  #   if(nrow(x) > 1){out = t.test(x$diff)$p.value}
  #   else{out = NA}
  #   return(out)
  # })
}
multiplot(p1, p2, cols=2)

#6600FF title = element_text(size=18),
#axis.text = element_text(size=18), 

# ###  Plot results per outcome
# p = list()
# for(threshold in thresholds){
#   A = subset(Rmeas, AUC.full > threshold | AUC.obs > threshold)
#   print(sprintf('%d of %d kept at threshold=%0.2f', nrow(A), nrow(Rmeas), threshold))
# 
# }


# # Plot AUCs on Rmeas, comparing meas vs. meas + pred
# p = ggplot(Rmeas, aes(x=AUC.obs, y=AUC.full)) +
#   geom_rect(data=r, aes(fill=category), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf,
#             alpha = 0.4) +
#   geom_point() +
#   geom_abline(intercept=0, slope=1) +
#   geom_hline(yintercept=0.5, col='DarkGrey', linetype=2) +
#   geom_vline(xintercept=0.5, col='DarkGrey', linetype=2) +
#   geom_text_repel(data=Rmeas, aes(x=AUC.obs, y=AUC.full, label=feature), color='black',
#                   size=2, segment.size = 0.3, segment.color='black', 
#                   box.padding = unit(0.5, "lines"),
#                   point.padding = unit(0.1, "lines")) +
#   facet_wrap(~outcome) +
#   scale_fill_brewer(palette='Set3') +
#   xlab('AUC using ONLY MEASURED signatures') +
#   ylab('AUC using MEASURED + IMPUTED signatures') +
#   ggtitle('Does addition of imputed signatures improve classification accuracy?') +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   xlim(c(0,1)) + ylim(c(0,1))
# print(p)
# #ggsave(plot=p, filename=paste0(resDir, '/Does_imputation_help.pdf'), width=15, height=12)

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


# # Subset to evaluations on full set of drugs
# Rfull = RemoveDfColumns(subset(R, eval == 'eval_all'), 'eval')
# Rfull = subset(Rfull, feature %in% c('max','mean') | subset == 'full')
# Rfull$feature_s = paste(Rfull$feature, Rfull$subset, sep='_')
# # Hm this one tries to order the features based on the median across outcomes
# out = lapply(split(Rfull, Rfull$outcome), function(dat){
#   ggplot(dat, aes(x=reorder(feature_s, -AUC, FUN=median), y=AUC, fill=subset)) + 
#     stat_summary(fun.y="median", geom="bar") + ggtitle(unique(dat$outcome)) +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab('') + guides(fill=FALSE)
#   })
# multiplot(plotlist=out, cols=5)