

numSigs = NumSigs(tensors$meas, 'cell')

#G = list()
#G$meas = ComputeGeneGeneCor(tensors$meas, cellSpecific=TRUE)
for(subset in c('cv', 'pred', 'merge')){
  #G[[subset]] = lapply(tensors[[subset]], function(tensor) ComputeGeneGeneCor(tensor, cellSpecific=TRUE))
  #GCP = lapply(G[[subset]], function(g) CorMatrixList(G$meas, g))
  #save(GCP, file=DataDir(sprintf('results/tsize/large/GCP_%s.RData', subset)))
  load(DataDir(sprintf('results/tsize/large/GCP_%s.RData', subset)))
  PlotGCP(GCP, subset=subset, numSigs, legend=FALSE, title=TRUE)
}

load(DataDir('results/tsize/large/GCP_cv.RData'))
PlotGCP(GCP, subset='cv', numSigs, legend=TRUE, title=FALSE, 
        file='preservation_of_ggc_cv_MAIN.tiff')

if(exists('OUTPUT')){
  OUTPUT$GCP_cv_mean = apply(as.matrix(as.data.frame(GCP)), 2, mean)
  OUTPUT$GCP_cv_sd = apply(as.matrix(as.data.frame(GCP)), 2, sd)
}