library(ggplot2)
library(dplyr)
library(plotrix)
library(R.matlab)
library(tidyr)

writeToFile = TRUE
addLegend = FALSE
useMethodColors = FALSE

main = list(mean='1D-Mean', mean2='2D-Mean', dnpp='DNPP', tensor='Tensor')

for(method in c('mean', 'mean2', 'dnpp', 'tensor')){
  
  true = as.vector(tensors$meas)
  est = as.vector(tensors$cv[[method]])
  cor = cor(true, est, use='pairwise.complete.obs', method='pearson')
  
    if(writeToFile){tiff(PlotDir(sprintf('expression_scatter_%s.tiff', method)))}
    par(mar=c(5,6,4,2))
    if(useMethodColors){
      color = GetMethodColor(method)
    }else{
      color = 'navajowhite4'
    }
    lim = c(-0.7,0.7)
    plot(c(-1,1), c(-1,1), col='grey', lwd=3, xlim=lim, ylim=lim, type='l', xlab='', ylab='', cex.axis=1.5)
    par(new=TRUE)
    smoothScatter(true, est, nbin=500, pch=20, col=NULL, transformation=function(x) x^.10,
                  colramp = colorRampPalette(c('white', color)), 
                  main=main[[method]], cex.main=2, 
                  xlab='true', ylab='predicted', cex.lab=2.0, cex.axis=1.5, xlim=lim, ylim=lim)
    usr = par('usr')
    text(usr[2], usr[3], bquote(r == .(round(cor, 2))), adj=c(1.3,-1), cex=2)
    if(addLegend){
      legend('topleft', c('perfect, slope=1',sprintf('actual, slope=%0.2f', slope)), 
             lwd=8, col=c('black', 'blue'), cex=1.8, bty='n')
    }
    if(writeToFile){dev.off()}
}


