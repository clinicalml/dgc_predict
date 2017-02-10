CheckMap = function(map, ndata){
  nprobe = length(unique(map$probe))
  stopifnot(nprobe <= ndata)
  stopifnot(length(unique(map$GeneID)) <= nprobe)
  stopifnot(nprobe == nrow(map))
  stopifnot(nprobe >= 0.1*ndata)
}

CheckDataLimits = function(dataframe){
  stopifnot(max(na.omit(dataframe)) < 20)
  stopifnot(min(na.omit(dataframe)) >= 0)
}

CheckNames = function(names){
  stopifnot(is.vector(names))
  stopifnot(is.character(names))
}