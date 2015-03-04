MedianNorm=function(Data){
	if(ncol(Data)==1)stop("Only 1 sample!")
  geomeans <- exp(rowMeans(log(Data)))
	apply(Data, 2, function(cnts) median((cnts/geomeans)[geomeans >  0]))
}
