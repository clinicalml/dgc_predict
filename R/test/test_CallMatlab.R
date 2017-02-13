TestStartMatlab(){
  # See TestMatlab below, this essentially tests whether Matlab has been started appropriately
}

TestMatlab = function(mat=get0('matlab')){
  working = FALSE
  if(!is.null(mat) && isOpen(mat)){
    evaluate(mat, 'a=2+2;')
    a = getVariable(mat, 'a')$a
    if(a == 4){
      working = TRUE
    }
  }
  return(working)
}

TestCrossValidateTensor(){}
TestCompleteTensor(){}
