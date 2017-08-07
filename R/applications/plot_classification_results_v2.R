
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

# Results where AUCs were both below this in the two comparisons are thrown out
threshold = 0.5 

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
Rmeas = FilterByAUC(Rmeas, threshold)
idx = Rmeas$category == 'ATC'
Rmeas$outcome[idx] = paste0('ATC ' , Rmeas$outcome[idx])

### First compute overall difference between the two groups
out = t.test(x=Rmeas$AUC.full, y=Rmeas$AUC.obs, paired=TRUE)


### Plot results per model
print(ggplot(Rmeas, aes(x=model, y=diff, group=model, fill=model)) + geom_boxplot() +
        ggtitle(sprintf('Deltas per model, threshold = %0.1f', threshold)) + stat_compare_means())
ggsave(PlotDir(sprintf('Deltas_per_model_threshold_%0.2f.pdf', threshold)))
p_model = lapply(split(Rmeas, Rmeas$model), function(x) t.test(x$AUC.full, x$AUC.obs, paired=TRUE)$p.value)


## Plot results per feature

# Add numSigs so that I can color boxplot accordingly
numSigs = c(MCF7=1505, VCAP=1368, PC3=1340, A375=1168, A549=1139,
            HA1E=1127, HT29=1022, HCC515=934, HEPG2=798, NPC=441)
Rmeas$num_sigs = numSigs[Rmeas$feature]

p_feature = sapply(split(Rmeas, Rmeas$feature), function(x) t.test(x$AUC.full, x$AUC.obs, paired=TRUE)$p.value)
adjP = p.adjust(p_feature, method='BH')

p1 = ggplot(melt(Rmeas), aes(x=reorder(feature, -diff, FUN=median), y=diff, group=feature, fill=num_sigs)) +
  geom_hline(yintercept = 0, color='grey', lwd=1) +
  geom_boxplot(alpha=0.9) +
  ggtitle(sprintf('Deltas per feature, threshold = %0.1f', threshold)) +
  theme_bw() + ylim(c(-0.5,0.5)) + scale_fill_gradientn(colours=brewer.pal(5,'YlOrRd')) +
  theme(text = element_text(size=18), axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size=20, hjust=0.5),
        legend.justification='bottom', legend.position='bottom') +
  labs(fill='# Signatures  \n  Measured\n', x='', y=expression(AUC [impute] - AUC [orig]),
       title='Improvement in AUC per Drug Signature Type')
 
# ggsave(PlotDir('DeltaAUC_per_Feature.svg'), height=7, width=8)



p2 =  ggplot(Rmeas, aes(x=reorder(outcome, -diff, FUN=median), y=diff, group=outcome, fill=category)) + ylim(c(-0.5,0.5)) +
  geom_hline(yintercept = 0, color='grey', lwd=1) + geom_boxplot(alpha=0.8) + scale_fill_manual(values=c('#CC9966','#009999')) +
  ggtitle(sprintf('Deltas per outcome, threshold = %0.1f', threshold)) + theme_bw() +
  labs(x='', y=expression(AUC [impute] - AUC [orig]), title='Improvement in AUC per Prediction Task',
       fill='Category') +
  theme(text = element_text(size=18), axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size=20, hjust=0.5),
        legend.justification='bottom', legend.position='bottom')

# ggsave(PlotDir('DeltaAUC_per_Outcome.svg'), height=7, width=8)

p_outcome = lapply(split(Rmeas, Rmeas$outcome), function(x) t.test(x$AUC.full, x$AUC.obs, paired=TRUE)$p.value)
adjp = p.adjust(p_outcome, method='BH')

# multiplot(p1, p2, cols=2)
