ComputeGCP = function(T_meas, T_pred_list, plot=FALSE, subset='cv', file=NA){
  G_meas = ComputeGeneGeneCor(T_meas, cellSpecific=TRUE)
  G_pred_list = lapply(T_pred_list, function(tensor) ComputeGeneGeneCor(tensor, cellSpecific=TRUE))
  GCP = lapply(G_pred_list, function(G) CorMatrixList(G_meas, G))
  
  if(!is.na(file)){
    save(GCP, file=file)
  }

  # Make barplots
  if(plot){
    PlotGCP(GCP, subset=subset)
  }
  
  return(GCP)
}

ComputeCSpPres = function(T_meas, T_pred_list, plot=FALSE, subset='cv', cs_true=NULL, debug=FALSE){
  
  if(is.null(cs_true)){
    cs_true = ComputeCellSpecificity(T_meas)$cs
  }
  
  if(debug){
    n=30
    cs_true = cs_true[1:n]
    T_pred_list = lapply(T_pred_list, function(tensor) tensor[1:n,,])
  }
  
  cs_pred = lapply(T_pred_list, function(tensor) ComputeCellSpecificity(tensor)$cs)

  if(plot){
    p = PlotCSP(cs_true, cs_pred, subset=subset)
  }else{
    p = NULL
  }
  
  CSp_cor = lapply(cs_pred, function(cs) cor(cs_true, cs, use='pairwise'))
  CSp_slope = lapply(cs_pred, function(cs){fit = lm(cs ~ 0 + cs_true); fit$coefficients})
  return(list(CSp_cor=CSp_cor, CSp_slope=CSp_slope, cs_pred, p=p))
}

ComputePCT_AllModes = function(T_meas, T_pred_list){
  
  list[nDrug,nGene,nCell] = dim(T_meas)

  dVec = 1:nDrug
  gVec = 1:nGene
  cVec = 1:nCell
  
  names(dVec) = dimnames(T_meas)[[1]]
  names(gVec) = dimnames(T_meas)[[2]]
  names(cVec) = dimnames(T_meas)[[3]]

  PCTd = list()
  PCTg = list()
  PCTc = list()
  
  for(method in names(T_pred_list)){
    T_pred = T_pred_list[[method]]
    PCTd[[method]] = sapply(dVec, function(i) cor(as.vector(T_meas[i,,]), as.vector(T_pred[i,,]), use='pairwise'))
    PCTg[[method]] = sapply(gVec, function(i) cor(as.vector(T_meas[,i,]), as.vector(T_pred[,i,]), use='pairwise'))
    PCTc[[method]] = sapply(cVec, function(i) cor(as.vector(T_meas[,,i]), as.vector(T_pred[,,i]), use='pairwise'))
  }
   
  return(list(PCTd=PCTd, PCTg=PCTg, PCTc=PCTc))
}

ComputePCT = function(T_meas, T_pred){
  list[x_meas, x_pred] = Tensor2Vec(T_meas, T_pred)
  return(cor(x_meas, x_pred))
}

ComputePCTPerSig = function(T1, T2, format='df'){
  stopifnot(identical(dim(T1), dim(T2)))
  nDrug = dim(T1)[1]
  nCell = dim(T1)[3]
  
  R = array(data=NA, dim=c(nDrug,nCell))
  P = R
  for(i in 1:nDrug){
    for(j in 1:nCell){
      if(!is.na(T1[i,1,j]) && !is.na(T2[i,1,j])){
        out = cor.test(T1[i,,j], T2[i,,j])
        R[i,j] = out$estimate
        P[i,j] = out$p.value
      }
    }
  }
  
  adjP = p.adjust(na.omit(as.vector(P)), method='BH')
  
  if(format == 'df'){
    R = na.omit(as.vector(R))
    P = na.omit(as.vector(P))
    out = data.frame(R=R, P=P, adjP=adjP)
  }else if(format == 'list'){
    out = list(R=R, P=P, adjP=adjP)
  }else{
    stop('unexpected value for format')
  }
  return(out)
}

Tensor2Vec = function(T1, T2){
  x1 = as.vector(T1)
  x2 = as.vector(T2)
  
  idx1 = na.action(na.omit(x1))
  idx2 = na.action(na.omit(x2))
  idxOmit = union(idx1, idx2)
  
  if(length(idxOmit) > 0){
    x1 = x1[-idxOmit]
    x2 = x2[-idxOmit]
  }
  
  return(list(x1=x1, x2=x2))
}
