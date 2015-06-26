
RankNorm=function(Data){
	if(ncol(Data)==1)stop("Only 1 sample!")

	RankData=apply(Data, 2, rank)
	SortData=apply(Data, 2, sort)
	SortMean=rowMeans(SortData)
	SortMean[SortMean==0]=1
	NormMatrix=sapply(1:ncol(Data), function(i)Data[,i]/(SortMean[RankData[,i]]))
	NormMatrix[NormMatrix==0]=1
	NormMatrix
	}

