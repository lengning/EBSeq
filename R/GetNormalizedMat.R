GetNormalizedMat<-function(Data, Sizes){
if(length(Sizes)!=length(Data) &  length(Sizes)!=ncol(Data))
	    stop("The number of library size factors is not the same as the number of samples!")
if(length(Sizes)==length(Data))Out=Data/Sizes
if(length(Sizes)==ncol(Data))Out=t(t(Data)/Sizes)
Out}
