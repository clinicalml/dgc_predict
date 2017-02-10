

numSigs = NumSigs(tensors$meas, 'cell')

# G = list()
# G$meas = ComputeGeneGeneCor(tensors$meas, cellSpecific=TRUE)
for(subset in c('cv', 'pred', 'merge')){
  #   #G[[subset]] = lapply(tensors[[subset]], function(tensor) ComputeGeneGeneCor(tensor, cellSpecific=TRUE))
  #   #GCP = lapply(G[[subset]], function(g) CorMatrixList(G$meas, g))
  #   #save(GCP, file=ResultsDir(sprintf('large/GCP_%s.RData', subset)))

  # For Supplementary
  load(ResultsDir(sprintf('large/GCP_%s.RData', subset)))
  PlotGCP(GCP, subset=subset, numSigs, legend=FALSE, title=TRUE)
}

# For Figure 5B
load(ResultsDir('large/GCP_cv.RData'))
PlotGCP(GCP, subset='cv', numSigs, legend=TRUE, title=FALSE, 
        file='preservation_of_ggc_cv_MAIN.tiff')
