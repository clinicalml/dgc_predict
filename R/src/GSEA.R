library(GSEABase)
library(HTSanalyzeR)
library(cellHTS2)
library(GO.db)
library(KEGG.db)
library(org.Hs.eg.db)
library(WriteXLS)
library(VennDiagram)
library(RCurl)
library(snow)
merge = base::merge

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
  
  #options(cluster = makeCluster(40, 'SOCK'))
  gsca = analyze(gsca, para=list(pValueCutoff=0.05, pAdjustMethod='BH', nPermutations=1000, 
                                 minGeneSetSize=minGeneSetSize, exponent=1), doGSOA=doGSOA, doGSEA=doGSEA)
  
  if(is(getOption('cluster'), 'cluster')) {
    stopCluster(getOption('cluster'))
    options(cluster=NULL)
  }
  
  gsca = AppendGSTerms(gsca)
  HTSanalyzeR::summarize(gsca)
  
  return(gsca)
}

AppendGSTerms = function(gsca, go = c('GO_BP', 'GO_MF', 'GO_CC'), 
                          kegg = c('PW_KEGG')){
  gsca = appendGSTerms(gsca, goGSCs=go, keggGSCs=kegg)
  return(gsca)
}

GetGSCAHits = function(gscaList, analysis='HyperGeo.results', gsc){
  out = RestrictColumns(GetGeneSet(gscaList[[1]], analysis, gsc), names(gscaList[1]))
  for(i in 2:length(gscaList)){
    tmp = RestrictColumns(GetGeneSet(gscaList[[i]], analysis, gsc), names(gscaList)[i])
    out = merge(out, tmp, by='Gene.Set.Term')
  }
  to.keep = which(sapply(1:nrow(out), function(i) length(which(
    out[i, 2:(length(gscaList)+1)] < 0.05))) > 0)
  return(out[to.keep, ])
}

RestrictColumns = function(gs, name){
  gs = gs[,c('Gene.Set.Term', 'Adjusted.Pvalue')]
  names(gs) = c('Gene.Set.Term', paste0('adj.p.', name))
  return(gs)
}

CatDF = function(dfList){
  if(length(dfList) == 0){
    DF = data.frame()
  }else if(length(dfList) == 1){
    DF = dfList[[1]]
  }
  if(length(dfList) > 1){
    DF =rbind(dfList[[1]], dfList[[2]])
  }
  if(length(dfList) > 2){
    for(i in 3:length(dfList)){
      DF = rbind(DF, dfList[[i]])
    }
  }
  return(DF)
}

GetTerm = function(gs){
  if(substr(gs, 1, 2) == 'GO'){
    term = eval(parse(text=paste0('Term(GOTERM$\'', gs,'\')')))
  }else if(substr(gs, 1, 2) == 'hs'){
    keggid2keggname = as.list(KEGGPATHID2NAME)
    term = keggid2keggname[substring(gs, 4)][[1]]
  }else{
    stop('error: unexpected gs')
  }
  return(term)
}

LoadGSCA = function(sigName, gs='msigdb'){
  cat(paste(sigName, '\n'))
  if(identical(gs, 'msigdb')){
    file = sprintf('results/gsea/cf_meta_v2/rdata/%s_gsca.RData', sigName)
  }else{
    file = sprintf('results/gsea/cf_meta_v2/rdata/%s_gsca_go_kegg.RData', sigName)
  }
  return(LoadGSCAFile(file, gs))
}

LoadGSCAFile = function(file, gs){
  load(file)
  if(!identical(gs, 'msigdb')){
    gsca = appendGSTerms(gsca, goGSCs=c('GO_BP', 'GO_MF', 'GO_CC'), keggGSCs=c('PW_KEGG'))
  }
  return(gsca)
}

WriteSigGeneSetsGoKegg = function(gsca, analysis, sigName){
  go_bp = OrderBy(GetGeneSet(gsca, analysis, 'GO_BP'), 'Adjusted.Pvalue')
  print(head(go_bp))
  go_mf = OrderBy(GetGeneSet(gsca, analysis, 'GO_MF'), 'Adjusted.Pvalue')
  go_cc = OrderBy(GetGeneSet(gsca, analysis, 'GO_CC'), 'Adjusted.Pvalue')
  pw_kegg = OrderBy(GetGeneSet(gsca, analysis, 'PW_KEGG'), 'Adjusted.Pvalue')
  file = sprintf('results/gsea/rescue/xls/%s_%s.xls', sigName, analysis)
  WriteXLS(c('go_bp', 'go_mf', 'go_cc', 'pw_kegg'), ExcelFileName=file, 
           SheetNames=as.vector(c('GO_BP','GO_MF','GO_CC','PW_KEGG')), AdjWidth=TRUE, BoldHeaderRow=TRUE)
}

GetGeneSet = function(gsca, analysis, oneGS){
  tableName = sprintf('gsca@result$%s$%s', analysis, oneGS)
  print(tableName)
  gsTable = eval(parse(text=tableName))
  if('FDR' %in% colnames(gsTable)){
    idx = which('FDR' == colnames(gsTable))
    gsTable = gsTable[,setdiff(1:ncol(gsTable), idx)]
  }
  return(gsTable)
}

OrderByFDR = function(table){
  stopifnot('FDR' %in% names(table))
  return(table[order(table$FDR),])
}

OrderBy = function(table, colName){
  table = as.data.frame(table)
  stopifnot(colName %in% names(table))
  expr = sprintf('table[order(table$%s),]', colName)
  orderedTable = eval(parse(text=expr))
  return(orderedTable)
}

TestOrderBy = function(){
  table = gsca1@result$GSEA.results$GO_MF
  orderedTable1 = OrderBy(table, 'FDR')
  orderedTable2 = OrderByFDR(table)
  stopifnot(identical(orderedTable1, orderedTable2))
}

GetCF2SigNames = function(){
  return(c('cf_meta_sam_fisher-05_eff-05_numstud-2_sorteff100',
           'cf_meta_sam_fisher-05_eff-05_numstud-2_sorteff200',
           'cf_meta_sam_fisher-05_eff-05_numstud-2_sorteff300',
           'cf_meta_fisherp1-05_eff-05_numstud-2'))
}

WriteGeneSetsToXLS = function(sigName, shortName, baseDir='results/gsea', metaSig=master14){
  load(paste0(baseDir, '/rdata/', sigName, '.RData'))
  gsca = appendGSTerms(gsca, goGSCs=c('GO_BP', 'GO_MF', 'GO_CC'), keggGSCs=c('PW_KEGG'))
  gscs = c('GO_BP', 'GO_MF', 'GO_CC', 'PW_KEGG')
  topGeneSets = getTopGeneSets(gsca, resultName='HyperGeo.results', gscs=gscs, allSig=TRUE)
  file = paste(baseDir, '/xls/', sigName, '_HyperGeo_GENE_SETS.xls', sep='')
  cat(paste0('\n', file, '\n'))
  
  colors = c('red', 'blue', 'green', 'yellow')
  df = list()
  genesInSigList = list()
  termList = list()
  for(gsc in names(topGeneSets)){
    for(gs in unlist(topGeneSets[gsc])){
      allGenesInSet = unlist(eval(parse(text=paste0('gsca@listOfGeneSetCollections$', 
                                                     gsc, '$\'', gs, '\''))))
      genesInSig = intersect(allGenesInSet, gsca@hits)
      genesInSigList[[length(genesInSigList)+1]] = genesInSig
      idx = match(genesInSig, metaSig$GeneID)
      namesInSig = metaSig[idx, 'Symbol']
      termList[[length(termList)+1]] = print(GetTerm(gs))
      
      df[[length(df) + 1]] = data.frame(GeneID=genesInSig, Name=namesInSig, Term=GetTerm(gs))
    }
    eval(parse(text=sprintf('%s = CatDF(df)', gsc)))
    
    file = PlotDir(paste(shortName, gsc, 'venn.tiff', sep='_'), subdir='gsea')
    
    # make venn diagrams
    if((length(genesInSigList) <= 5) & (length(genesInSigList) > 0)){
      names(genesInSigList) = termList 
      venn.diagram(genesInSigList, file=file,
                   main=paste(shortName, gsc, 'hits', sep='_'), 
                   fill=rep_len(colors, length.out=length(genesInSigList))) 
    }
    df = list()
    genesInSigList = list()
    termList = list()
  }
  
  #WriteXLS(names(topGeneSets), ExcelFileName=file, names(topGeneSets), AdjWidth=FALSE, BoldHeaderRow=TRUE)  
}

GetGSEAScores = function(gscaList){
  out = list()
  for(i in 1:length(gscaList)){
    gsca = names(gscaList)[i]
    R = gscaList[[i]]@result$GSEA.results

    if(length(R) > 1){
      out[[gsca]] = list()
      for(geneset in names(R)){
        out[[gsca]][[geneset]] = R[[geneset]]
      }
    }else{
      out[[gsca]] = R[[1]]
    }
  }
  return(out)
}

GetGSEAResultsVector = function(singleGSEAScores, colName, 
                                 genesets=c('GO_MF', 'GO_BP', 'GO_CC', 'PW_KEGG')){
  if(genesets == 'all'){
    genesets = names(singleGSEAScores)
  }
  stopifnot(genesets %in% names(singleGSEAScores))
  x = c()
  for(geneset in genesets){
    x = c(x,singleGSEAScores[[geneset]][,colName])
  }
  return(x)
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

# try pearson, spearman correlation, Fisher's exact test, cosine distance, edit distance
GSEACompare = function(x, y, metric='pearson', thresh=0.05){
  x = na.omit(x)
  y = na.omit(y)
  
  toKeep = intersect(names(x), names(y))
  x = x[toKeep]
  y = y[toKeep]
  
  if(length(x) == 0 || length(y) == 0){
    warning('x or y is all NA')
    out = NA
  }else{
    if(metric == 'pearson'){
      out = cor(x,y)
    }else if(metric == 'spearman'){
      out = cor(x,y,method='spearman')
    }else if(metric == 'fisher'){
      A = names(x)[which(x < thresh)]
      B = names(y)[which(y < thresh)]
      out = -log10(FisherExact(A=A, B=B, n_universe=length(x), printFlag = FALSE)$p)
    }else if(metric == 'hinge'){
      xx = sapply(-log10(x), function(i) min(i,3))
      yy = sapply(-log10(y), function(i) min(i,3))
      out = cor(xx,yy)
    }else if(metric == 'cos'){
      out = sum(x*y) / (Norm2(x)*Norm2(y))
    }else{
      warning('did not recognize metric')
    }
  }
  return(out)
}

SplitB = function(B){
  control = as.vector(B[1:4,])
  matched = c(B[5,1], B[6,2], B[7,3], B[8,4], B[9,5], B[10,6])
  return(list(control=control, matched=matched))
}

ExtractGSOAResults = function(gsca, thresh=1){
  out = list()
  for(type in c('GO_MF',  'GO_BP', 'GO_CC', 'PW_KEGG')){
    tmp = gsca@result$HyperGeo.results[[type]]
    A = tmp[tmp$Pvalue<thresh,]
    A$geneSet = type
    A = A[order(A$Pvalue),]
    out[[type]] = A
  }
  return(Reduce('rbind', out))
}
