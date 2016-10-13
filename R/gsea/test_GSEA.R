

TestRunGsea <- function(){
  warning('needs to be further tested')
  
  set.seed(123)
  genes = as.character(GetLmGenes(type='entrez'))[1:400]
  x = rnorm(length(genes))
  names(x) = genes
  hits = sample(genes, 10)
  
  # check that ordering of x doesn't make a difference when running GSEA
  set.seed(123)
  test1 = RunGsea(profile=x, dz_genes_up=hits, minGeneSetSize=20, doGSOA=FALSE, doGSEA=TRUE)
  set.seed(123)
  test2 = RunGsea(profile=sort(x), dz_genes_up=hits, minGeneSetSize=20, doGSOA=FALSE, doGSEA=TRUE)
  stopifnot(identical(test1@result$GSEA.results, test2@result$GSEA.results))
  
  
  # check that sorting x doesn't make a difference when running GSOA
  set.seed(123)
  test1 = RunGsea(profile=x, dz_genes_up=hits, minGeneSetSize=20, doGSOA=TRUE, doGSEA=FALSE)
  set.seed(123)
  test2 = RunGsea(profile=sort(x), dz_genes_up=hits, minGeneSetSize=20, doGSOA=TRUE, doGSEA=FALSE)
  stopifnot(identical(test1@result$HyperGeo.results, test2@result$HyperGeo.results))
}