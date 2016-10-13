# This code evaluates AUCs as a function of various ways to call DEGs across the
# tensor.  DEGs are always called based on some threshold (so that abs(expr) >
# thresh is called differentially expressed), however this threshold can be
# defined differently over various subsets of the tensor, and some preprocessing
# can be done.  I’m playing with the following knobs here:
# 
# 1)	‘method’ as either tensor, sig, or gene.  Tensor works with all values of 
# the tensor at once, whereas sig looks at individual signatures (i.e. a single 
# drug/cell combination), and gene looks at all values for that particular gene.
# 2)	‘ percDEG’ specifies the desired percentage of genes to consider 
# significantly differentially expressed.  I looked a range between 0.5-10% 
# 3)	‘symmetric’ specifies whether I should call an equal number of up-regulated
# and down-regulated genes, or whether I should use a fixed threshold, and hence
# allow for an asymmetric number of up vs. down genes.
# 4)	‘normGene’ determines whether both input and true tensor should be 
# preprocessed by centering the distribution of values for each gene to a standard
# normal (i.e. converting gene values to z-scores.)
# 
# debug = FALSE
# 
# # load data
# Tmeas = list(nonorm=GetDataTensor())
# Tcv = list(nonorm=GetAllTensors(meas=FALSE, cv=TRUE, comp=FALSE, merge=FALSE, pred=FALSE)$cv)
# if(debug){
#   d = 30
#   Tcv$nonorm$mean  = Tcv$nonorm$mean[1:d,,]
#   Tcv$nonorm$mean2  = Tcv$nonorm$mean2[1:d,,]
#   Tcv$nonorm$tensor = Tcv$nonorm$tensor[1:d,,]
#   Tmeas$nonorm  = Tmeas$nonorm[1:d,,]
# }
# 
# # Tmeas$norm = TensorZ(Tmeas$nonorm)$Z
# # Tcv$norm$mean = TensorZ(Tcv$nonorm$mean)$Z
# # Tcv$norm$mean2 = TensorZ(Tcv$nonorm$mean2)$Z
# # Tcv$norm$tensor = TensorZ(Tcv$nonorm$tensor)$Z
# 
# slices = c('tensor', 'sig')#, 'gene')
# percDEG = c(0.5, 1, 2, 5, 10)
# syms = c(FALSE, TRUE)
# normStr = list(nonorm=FALSE) #norm=TRUE, 
# 
# dim = c(length(slices), length(percDEG), length(syms), length(normStr))
# 
# dimnames = list(slice =slices, 
#                 percDEG = as.character(percDEG),
#                 symmetric = as.character(syms),
#                 normGene = names(normStr))
# 
# tmp = array(data=NA, dim=dim, dimnames=dimnames)
# auc = list(mean=tmp, mean2=tmp, tensor=tmp)
# 
# AUC = data.frame(slice='', percDEG=0, symmetric=FALSE, normGene='', tcMethod='', value=0)
# 
# for(slice in slices){
#   print(sprintf('slice=%s', slice))
#   for(perc in percDEG){
#     print(sprintf(' perc=%.2f', perc))
#     for(sym in syms){
#       print(sprintf('  sym=%s', sym))
#       for(norm in names(normStr)){
#         print(sprintf('   normGene=%s', norm))
#         
#         D = TensorDEG(Tmeas[[norm]], normGene=normStr[[norm]], method=slice, percDEG=perc, symmetric=sym)$D
#         labels = as.vector(abs(D))
#         
#         for(tcMethod in c('mean', 'mean2', 'tensor')){
#           val = ComputeAUC(est=as.vector(Tcv[[norm]][[tcMethod]]), labels=labels, computeROC=FALSE)
#           AUC = rbind(AUC, c(slice, perc, sym, norm, tcMethod, val))
#         }
#       }
#     }
#   }
# }
# 
# AUC = AUC[2:nrow(AUC),]
# class(AUC$symmetric) = 'logical'
# class(AUC$percDEG) = 'numeric'
# class(AUC$value) = 'numeric'
# save(AUC, file=ResultsDir('choose_tensorDEG_method.RData'))

# # plot
# for(norm in c('norm', 'nonorm')){
#   for(sym in c(TRUE, FALSE)){
#     auc = subset(AUC, symmetric==sym & normGene==norm)
#     p = ggplot(auc, aes(x=percDEG, y=value, group=interaction(slice, tcMethod), color=tcMethod, linetype=slice)) +
#       geom_line(size=2) + 
#       xlab("percent DEG") +
#       ylab("AUC") +
#       ggtitle(sprintf('Sym=%s, NormGene=%s', sym, norm)) + 
#       guides(lty=guide_legend(keywidth=3), color=guide_legend(keywidth=3)) +
#       theme_classic() +
#       theme(text=element_text(size=26),
#             legend.title = element_blank()) +
#       scale_color_manual(values=GetMethodColors()) +
#       ylim(c(0.5,1))
#     print(p)
#   }
# }

#AUC$tcMethod = revalue(AUC$tcMethod, c('mean'='1D-Mean', 'mean2'='2D-Mean', 'tensor'='Tensor'))
#AUC$slice[AUC$slice == 'tensor'] = 'uniform'
#AUC$slice[AUC$slice == 'sig'] = 'per sig'

# compare sig vs. tensor
auc = subset(AUC, symmetric==FALSE)
p1 = ggplot(auc, aes(x=percDEG, y=value, 
                    group=interaction(slice, tcMethod),
                    color=tcMethod, linetype=slice)) +
  geom_line(size=1.5) + 
  xlab('% genes called as DEG') +
  ylab('AUC') +
  guides(lty=guide_legend(keywidth=2), color=guide_legend(keywidth=2)) +
  theme_classic() +
  theme(text=element_text(size=22), legend.position='bottom') +
  scale_color_manual(values=GetMethodColors(longName=TRUE)) +
  labs(linetype='threshold choice', color='method') +
  ylim(c(0.6,0.9))
tiff(PlotDir('DEG_slice.tiff'), width=550)
print(p1)
dev.off()

# compare symmetric vs. asymetric
auc = subset(AUC, slice=='per sig')
p2 = ggplot(auc, aes(x=percDEG, y=value, 
                    group=interaction(symmetric, tcMethod),
                    color=tcMethod, linetype=symmetric)) +
  geom_line(size=1.5) + 
  xlab('% genes called as DEG') +
  ylab('AUC') +
  guides(lty=guide_legend(keywidth=2), color=guide_legend(keywidth=2)) +
  theme_classic() +
  theme(text=element_text(size=22), legend.position='bottom') +
  scale_color_manual(values=GetMethodColors(longName=TRUE)) +
  labs(linetype='symmetric', color='method') +
  ylim(c(0.6,0.9))
print(p2)

tiff(PlotDir('DEG_symmetric.tiff'), width=550)
print(p2)
dev.off()
