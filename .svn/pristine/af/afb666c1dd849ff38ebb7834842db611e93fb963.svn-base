GetMultiFC=function(EBMultiOut,SmallNum=.01){
	if(!"PPpattern"%in%names(EBMultiOut))stop("The input doesn't seem like an output from EBMultiTest")

	NumNgGroup=length(EBMultiOut$DataList)
	OutNames=rownames(EBMultiOut$PPMat)
	NumCondition=length(EBMultiOut$SPMean[[1]])
	ConditionNames=colnames(EBMultiOut$AllParti)
	CondMeans=sapply(1:NumCondition,
							function(i){
								if (NumNgGroup==1)
									out=EBMultiOut$SPMean[[1]][[i]][OutNames]
								if(NumNgGroup>1)
								out=unlist(sapply(1:NumNgGroup,
							function(j)EBMultiOut$SPMean[[j]][[i]]))[OutNames]	
							out}
								)
	colnames(CondMeans)=ConditionNames
	CondMeansPlus=CondMeans+SmallNum

	GeneRealMean=rowMeans(CondMeans)
	GeneR=unlist(EBMultiOut$RList)
  GeneR[GeneR<=0 | is.na(GeneR)]=GeneRealMean[GeneR<=0 | is.na(GeneR)]*.99/.01

  GeneAlpha=EBMultiOut[[1]][nrow(EBMultiOut[[1]]),]
  GeneBeta=unlist(sapply(1:length(EBMultiOut$DataList),
												 function(i)rep(EBMultiOut[[2]][nrow(EBMultiOut[[1]]),i],
																				nrow(EBMultiOut$DataList[[i]]))))
  GeneBeta=as.vector(GeneBeta)


	FCMat=PostFCMat=matrix(0,ncol=choose(NumCondition,2),nrow=length(OutNames))
	rownames(FCMat)=rownames(PostFCMat)=OutNames
	k=1
	ColNames=rep(NA,choose(NumCondition,2))
	for(i in 1:(NumCondition-1)){
		for(j in (i+1):NumCondition)
		{
		ColNames[k]=paste(ConditionNames[i],"Over",ConditionNames[j],sep="")
		FCMat[,k]=CondMeansPlus[,i]/CondMeansPlus[,j]


		nC1=sum(EBMultiOut$ConditionOrder==ConditionNames[i])
		nC2=sum(EBMultiOut$ConditionOrder==ConditionNames[j])
		GenePostAlphaC1=GeneAlpha+nC1*GeneR
		GenePostAlphaC2=GeneAlpha+nC2*GeneR
		GenePostBetaC1=GeneBeta+nC1*CondMeans[,i]
		GenePostBetaC2=GeneBeta+nC2*CondMeans[,j]
		GenePostQC1=GenePostAlphaC1/(GenePostAlphaC1+GenePostBetaC1)
		GenePostQC2=GenePostAlphaC2/(GenePostAlphaC2+GenePostBetaC2)

		GenePostFC=((1-GenePostQC1)/(1-GenePostQC2))*(GenePostQC2/GenePostQC1)
		PostFCMat[,k]= GenePostFC

		k=k+1
		}
	}
	colnames(FCMat)=colnames(PostFCMat)=ColNames
	Log2FCMat=log2(FCMat)
	Log2PostFCMat=log2(PostFCMat)
	Out=list(FCMat=FCMat,Log2FCMat=Log2FCMat,
					 PostFCMat=PostFCMat, Log2PostFCMat=Log2PostFCMat,
					 CondMeans=CondMeans, 
					 ConditionOrder=EBMultiOut$ConditionOrder)
}

