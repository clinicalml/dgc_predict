
RunGsea = function(profile, dz_genes_up, dz_genes_down=NULL, minGeneSetSize=15, 
                    doGSOA=FALSE, doGSEA=TRUE){
  
  GO_MF = GOGeneSets(species='Hs', ontologies=c('MF'))
  GO_BP = GOGeneSets(species='Hs', ontologies=c('BP'))
  GO_CC = GOGeneSets(species='Hs', ontologies=c('CC'))
  PW_KEGG = KeggGeneSets(species='Hs')
  ListGSC = list(GO_MF=GO_MF, GO_BP=GO_BP, GO_CC=GO_CC, PW_KEGG=PW_KEGG)
  
  cat('initializing gsca object...\n')
  gsca = new('GSCA', listOfGeneSetCollections=ListGSC, geneList=profile, hits=as.character(c(dz_genes_up, dz_genes_down)))
  
  cat('preprocessing...\n')
  gsca = preprocess(gsca, species='Hs', initialIDs='Entrez.gene', keepMultipleMappings=TRUE, 
                    duplicateRemoverMethod='max', orderAbsValue=FALSE)
  
  cat('analyzing...\n')

  gsca = analyze(gsca, para=list(pValueCutoff=0.05, pAdjustMethod='BH', nPermutations=1000, 
                                 minGeneSetSize=minGeneSetSize, exponent=1), doGSOA=doGSOA, doGSEA=doGSEA)
  gsca = AppendGSTerms(gsca)
  HTSanalyzeR::summarize(gsca)
  return(gsca)
}

AppendGSTerms = function(gsca, go = c('GO_BP', 'GO_MF', 'GO_CC'), kegg = c('PW_KEGG')){
  return(appendGSTerms(gsca, goGSCs=go, keggGSCs=kegg))
}

GetGSEAResultsV2 = function(gsca, colname='Pvalue', thresh=0.01, analysis='GSOA', merge=FALSE){
  if(analysis == 'GSOA'){
    x = gsca@result$HyperGeo.results
  }else if(analysis == 'GSEA'){
    x = gsca@result$GSEA.results
  }else{
    stop('unexpected value for analysis')
  }
  
  out = list()
  for(n in names(x)){
    y = x[[n]][,colname]
    out[[n]] = x[[n]][y <= thresh,]
  }
  
  if(merge){
    geneSets = names(which(sapply(out, nrow) > 0))
    out = lapply(geneSets, function(name){out[[name]]$Gene.Set=name; return(out[[name]])})
    out = Reduce('rbind', out)
  }
  return(out)
}