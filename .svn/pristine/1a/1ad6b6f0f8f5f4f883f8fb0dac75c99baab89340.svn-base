GetPatterns<-function(Conditions){
  if(!is.factor(Conditions))Conditions=as.factor(Conditions)
	NumCond=nlevels(Conditions)
	if(NumCond<3)stop("Less than 3 conditions!")
	CondLevels=levels(Conditions)
    AllPartiList=sapply(1:NumCond,function(i)nkpartitions(NumCond,i))
    AllParti=do.call(rbind,AllPartiList)
	colnames(AllParti)=CondLevels
	rownames(AllParti)=paste("Pattern",1:nrow(AllParti),sep="")
	AllParti

}
