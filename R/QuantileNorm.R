
QuantileNorm=function(Data, Quantile){
	if(ncol(Data)==1)stop("Only 1 sample!")

	QtilePt=apply(Data, 2, function(i)quantile(i, Quantile))
#	Size= QtilePt * prod(QtilePt) ^ (-1/ncol(Data))
	Size=10^(log10(QtilePt)-sum(log10(QtilePt))*(1/ncol(Data)) )
	Size
	}

