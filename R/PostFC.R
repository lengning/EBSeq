PostFC=function(EBoutput, SmallNum=.01) {
	if(!"C1Mean"%in%names(EBoutput))
		stop("The input doesn't seem like an output from EBTest")
	GeneRealMeanC1=unlist(EBoutput$C1Mean)
	GeneRealMeanC2=unlist(EBoutput$C2Mean)
	GeneRealMeanC1Plus=GeneRealMeanC1+SmallNum
	GeneRealMeanC2Plus=GeneRealMeanC2+SmallNum
	GeneRealMean=(GeneRealMeanC1+GeneRealMeanC2)/2

	GeneRealFC=GeneRealMeanC1Plus/GeneRealMeanC2Plus

	GeneR=unlist(EBoutput$RList)
	GeneR[GeneR<=0 | is.na(GeneR)]=GeneRealMean[GeneR<=0 | is.na(GeneR)]*.99/.01

	GeneAlpha=EBoutput[[1]][nrow(EBoutput[[1]]),]
	GeneBeta=unlist(sapply(1:length(EBoutput$C1Mean),function(i)rep(EBoutput[[2]][nrow(EBoutput[[1]]),i],length(EBoutput$C1Mean[[i]]))))
	GeneBeta=as.vector(GeneBeta)
	  
		# Post alpha P_a_C1= alpha + r_C1 * n_C1
	  # Post beta P_b_C1= beta + Mean_C1 * n_C1
	  # P_q_C1= P_a_C1/ (P_a_C1 + P_b_C1)
	  # Post FC = ((1-P_q_C1)/P_q_c1) /( (1-P_q_c2)/P_q_c2)

	nC1=sum(EBoutput$Conditions==levels(EBoutput$Conditions)[1])
	nC2=sum(EBoutput$Conditions==levels(EBoutput$Conditions)[2])
	GenePostAlphaC1=GeneAlpha+nC1*GeneR
	GenePostAlphaC2=GeneAlpha+nC2*GeneR
	GenePostBetaC1=GeneBeta+nC1*GeneRealMeanC1
	GenePostBetaC2=GeneBeta+nC2*GeneRealMeanC2
	GenePostQC1=GenePostAlphaC1/(GenePostAlphaC1+GenePostBetaC1)
	GenePostQC2=GenePostAlphaC2/(GenePostAlphaC2+GenePostBetaC2)

	GenePostFC=((1-GenePostQC1)/(1-GenePostQC2))*(GenePostQC2/GenePostQC1)
	Out=list(PostFC=GenePostFC[rownames(EBoutput$PPMat)], RealFC=GeneRealFC[rownames(EBoutput$PPMat)],
					 Direction=paste(EBoutput$ConditionOrder[[1]],"Over", EBoutput$ConditionOrder[[2]])
					 )

}
