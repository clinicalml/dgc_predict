TestRunGsea = function(){
  lmGenes = GetLmGenes('entrez')
  
  GO_CC = GOGeneSets(species='Hs', ontologies=c('CC'))
  PW_KEGG = KeggGeneSets(species='Hs')
  ListGSC = list(GO_CC=GO_CC, PW_KEGG=PW_KEGG)
  
  geneSet = intersect(union(GO_CC[[480]],PW_KEGG[[52]]), lmGenes)[-1]
  x = rnorm(978)
  names(x) = lmGenes
  
  capture.output(gsea <- RunGsea(profile=x, dz_genes_up=geneSet, minGeneSetSize=13, 
                                 doGSOA=TRUE, doGSEA=FALSE, ListGSC=ListGSC), file='/dev/null')
  A = gsea@result$HyperGeo.results$GO_CC
  stopifnot(A[order(A$Adjusted.Pvalue),'Gene.Set.Term'][1] == 'PML body')
  
  res1 = GetGSEAResultsV2(gsea, merge = FALSE)
  res2 = GetGSEAResultsV2(gsea, merge = TRUE)
  stopifnot(nrow(res2) == sum(sapply(res1, nrow)))
  
  stopifnot(res1$GO_CC$Gene.Set.Term[1] == 'PML body')
  stopifnot(res1$PW_KEGG$Gene.Set.Term[1] == 'Leishmaniasis')
}

TestAppendGSTerms = function(){} # Tested by TestRunGSEA 
TestGetGSEAResultsV2 = function(){} # Tested by TestRunGSEA 
TestGetListGSC = function(){}
