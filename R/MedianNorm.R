MedianNorm <- function(Data, alternative=FALSE){
	if(ncol(Data)==1)stop("Only 1 sample!")
  if(!alternative){
		geomeans <- exp(rowMeans(log(Data)))
		out <- apply(Data, 2, function(cnts) median((cnts/geomeans)[geomeans >  0]))}
	if(alternative){
		DataMatO <- Data
		N <- ncol(DataMatO)
		DataList0 <- sapply(1:N,function(i)DataMatO[,i]/DataMatO,simplify=F)
		DataEachMed0 <- sapply(1:N,function(i)apply(DataList0[[i]],2,function(j)median(j[which(j>0
		& j<Inf)])))
		DataColgeo <- sapply(1:N,function(i)exp(mean(log(DataEachMed0[-i,i]))))
		out <- DataColgeo
		}

	out
}
