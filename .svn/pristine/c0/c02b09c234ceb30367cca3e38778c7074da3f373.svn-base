EBTest <-
function(Data,NgVector=NULL,Conditions, sizeFactors, maxround, Pool=F, NumBin=1000,ApproxVal=10^-10, Alpha=NULL, Beta=NULL,PInput=NULL,RInput=NULL,PoolLower=.25, PoolUpper=.75,Print=T, Qtrm=.75,QtrmCut=10)
{
	if(!is.factor(Conditions))Conditions=as.factor(Conditions)
	if(is.null(rownames(Data)))stop("Please add gene/isoform names to the data matrix")

	if(!is.matrix(Data))stop("The input Data is not a matrix")
	if(length(Conditions)!=ncol(Data))stop("The number of conditions is not the same as the number of samples! ")
	if(nlevels(Conditions)>2)stop("More than 2 conditions! Please use EBMultiTest() function")
	if(nlevels(Conditions)<2)stop("Less than 2 conditions - Please check your input")
	if(length(sizeFactors)!=length(Data) &  length(sizeFactors)!=ncol(Data))
		stop("The number of library size factors is not the same as the number of samples!")		
	
	Conditions=as.factor(Conditions)
	Vect5End=Vect3End=CI=CIthre=tau=NULL
	Dataraw=Data

	#Normalized
	DataNorm=GetNormalizedMat(Data, sizeFactors)
	Levels=levels(as.factor(Conditions))

	# Dixon Statistics
#	library(outliers)
# normalized matrix for each condition
#	matC=sapply(1:length(Levels),function(i)DataNorm[,which(Conditions==Levels[i])])
# run dixon test for each isoform within condition
#	DixonP=sapply(1:length(matC),function(j)
#	apply(DataNorm,1,function(i){
#		  if(mean(i)==0)out=NA
#				  else out=dixon.test(i)$p.value
#						    out}))


	QuantileFor0=apply(DataNorm,1,function(i)quantile(i,Qtrm))
	AllZeroNames=which(QuantileFor0<=QtrmCut)
	NotAllZeroNames=which(QuantileFor0>QtrmCut)
	if(length(AllZeroNames)>0 & Print==T)
					    cat(paste0("Removing transcripts with ",Qtrm*100,
							    " th quantile < = ",QtrmCut," \n",
									length(NotAllZeroNames)," transcripts will be tested\n"))
	if(length(NotAllZeroNames)==0)stop("0 transcript passed")
	Data=Data[NotAllZeroNames,]
	if(!is.null(NgVector))NgVector=NgVector[NotAllZeroNames]
	if(length(sizeFactors)==length(Data))sizeFactors=sizeFactors[NotAllZeroNames,]
	if(is.null(NgVector))NgVector=rep(1,nrow(Data))

	#Rename Them
	IsoNamesIn=rownames(Data)
	Names=paste("I",c(1:dim(Data)[1]),sep="")
	names(IsoNamesIn)=Names
	rownames(Data)=paste("I",c(1:dim(Data)[1]),sep="")
	names(NgVector)=paste("I",c(1:dim(Data)[1]),sep="")
	

	if(length(sizeFactors)==length(Data)){
		rownames(sizeFactors)=rownames(Data)
		colnames(sizeFactors)=Conditions
	}
	
	NumOfNg=nlevels(as.factor(NgVector))
	NameList=sapply(1:NumOfNg,function(i)Names[NgVector==i],simplify=F)
	names(NameList)=paste("Ng",c(1:NumOfNg),sep="")
	NotNone=NULL
	for (i in 1:NumOfNg) {
		if (length(NameList[[i]])!=0) 
			NotNone=c(NotNone,names(NameList)[i])
		}
	NameList=NameList[NotNone]
		
	NoneZeroLength=length(NameList)
	DataList=vector("list",NoneZeroLength)
	DataList=sapply(1:NoneZeroLength , function(i) Data[NameList[[i]],],simplify=F)
	names(DataList)=names(NameList)
    
	NumEachGroup=sapply(1:NoneZeroLength , function(i)dim(DataList[[i]])[1])
	# Unlist 
	DataList.unlist=do.call(rbind, DataList)

	# Divide by SampleSize factor
	
	if(length(sizeFactors)==ncol(Data))
	DataList.unlist.dvd=t(t( DataList.unlist)/sizeFactors)
	
	if(length(sizeFactors)==length(Data))
	DataList.unlist.dvd=DataList.unlist/sizeFactors
	
	MeanList=rowMeans(DataList.unlist.dvd)

###############
# Input R
###############
	if (!is.null(RInput)){

	RNoZero=RInput[NotAllZeroNames]
	names(RNoZero)=rownames(Data)
	RNoZero.order=RNoZero[rownames(DataList.unlist)]
	if(length(sizeFactors)==ncol(Data)){
		RMat= outer(RNoZero.order, sizeFactors)
	}
	if(length(sizeFactors)==length(Data)){
		RMat= RNoZero.order* sizeFactors
	}
		
	DataListSP=vector("list",nlevels(Conditions))
	RMatSP=vector("list",nlevels(Conditions))

	for (lv in 1:nlevels(Conditions)){
		DataListSP[[lv]]= matrix(DataList.unlist[,Conditions==levels(Conditions)[lv]],nrow=dim(DataList.unlist)[1])
		rownames(DataListSP[[lv]])=rownames(DataList.unlist)
		RMatSP[[lv]]= matrix(RMat[,Conditions==levels(Conditions)[lv]],nrow=dim(RMat)[1])
		rownames(RMatSP[[lv]])=rownames(RMat)
	
		}

	F0Log=f0(Input=DataList.unlist, AlphaIn=Alpha, BetaIn=Beta, 
			 EmpiricalR=RMat, NumOfGroups=NumEachGroup, log=T)
	F1Log=f1(Input1=DataListSP[[1]], Input2=DataListSP[[2]], 
			 AlphaIn=Alpha, BetaIn=Beta, EmpiricalRSP1=RMatSP[[1]], 
			 EmpiricalRSP2=RMatSP[[2]], NumOfGroup=NumEachGroup, log=T)
	F0LogMdf=F0Log+600
	F1LogMdf=F1Log+600
	F0Mdf=exp(F0LogMdf)
	F1Mdf=exp(F1LogMdf)
	if(!is.null(PInput)){
		z.list=PInput*F1Mdf/(PInput*F1Mdf+(1-PInput)*F0Mdf)
		PIn=PInput
	}
	if(is.null(PInput)){
		PIn=.5
		PInput=rep(NULL,maxround)
		for(i in 1:maxround){
			z.list=PIn*F1Mdf/(PIn*F1Mdf+(1-PIn)*F0Mdf)
			zNaNName=names(z.list)[is.na(z.list)]
			zGood=which(!is.na(z.list))
			PIn=sum(z.list[zGood])/length(z.list[zGood])
			PInput[i]=PIn
		}
	
	zNaNName=names(z.list)[is.na(z.list)]
	if(length(zNaNName)!=0){
	PNotIn=rep(1-ApproxVal,length(zNaNName))
	MeanList.NotIn=MeanList[zNaNName]
	R.NotIn.raw=MeanList.NotIn*PNotIn/(1-PNotIn)
	if(length(sizeFactors)==ncol(Data))
			R.NotIn=outer(R.NotIn.raw,sizeFactors)
	if(length(sizeFactors)==length(Data))
			R.NotIn=R.NotIn.raw*sizeFactors[zNaNName,]
	R.NotIn1=matrix(R.NotIn[,Conditions==levels(Conditions)[1]],nrow=nrow(R.NotIn))
    R.NotIn2=matrix(R.NotIn[,Conditions==levels(Conditions)[2]],nrow=nrow(R.NotIn))
	NumOfEachGroupNA=sapply(1:NoneZeroLength, function(i)sum(zNaNName%in%rownames(DataList[[i]])))
	F0LogNA=f0(matrix(DataList.unlist[zNaNName,],ncol=ncol(DataList.unlist)),  Alpha, Beta, R.NotIn, NumOfEachGroupNA, log=T)
    F1LogNA=f1(matrix(DataListSP[[1]][zNaNName,],ncol=ncol(DataListSP[[1]])), 
			   matrix(DataListSP[[2]][zNaNName,],ncol=ncol(DataListSP[[2]])),
			   Alpha, Beta, R.NotIn1,R.NotIn2, NumOfEachGroupNA, log=T)
	F0LogMdfNA=F0LogNA+600
		F1LogMdfNA=F1LogNA+600
		F0MdfNA=exp(F0LogMdfNA)
		F1MdfNA=exp(F1LogMdfNA)
		z.list.NotIn=PIn*F1MdfNA/(PIn*F1MdfNA+(1-PIn)*F0MdfNA)
	z.list[zNaNName]=z.list.NotIn
	F0Log[zNaNName]=F0LogNA
	F1Log[zNaNName]=F1LogNA
	}
	}
	RealName.Z.output=z.list
	RealName.F0=F0Log
	RealName.F1=F1Log
	names(RealName.Z.output)=IsoNamesIn
	names(RealName.F0)=IsoNamesIn
	names(RealName.F1)=IsoNamesIn


	output=list(Alpha=Alpha,Beta=Beta,P=PInput, Z=RealName.Z.output,
			 PPDE=RealName.Z.output,f0=RealName.F0, f1=RealName.F1)
	return(output)
	
	}


	# Get FC and VarPool for pooling - Only works on 2 conditions
	if(ncol(Data)==2){
	DataforPoolSP.dvd1=matrix(DataList.unlist.dvd[,Conditions==levels(Conditions)[1]],nrow=dim(DataList.unlist)[1])	
	DataforPoolSP.dvd2=matrix(DataList.unlist.dvd[,Conditions==levels(Conditions)[2]],nrow=dim(DataList.unlist)[1])
	MeanforPoolSP.dvd1=rowMeans(DataforPoolSP.dvd1)
	MeanforPoolSP.dvd2=rowMeans(DataforPoolSP.dvd2)
	FCforPool=MeanforPoolSP.dvd1/MeanforPoolSP.dvd2
	names(FCforPool)=rownames(Data)
	FC_Use=which(FCforPool>=quantile(FCforPool[!is.na(FCforPool)],PoolLower) & 
								  FCforPool<=quantile(FCforPool[!is.na(FCforPool)],PoolUpper))
	
	Var_FC_Use=apply( DataList.unlist.dvd[FC_Use,],1,var )
	Mean_FC_Use=(MeanforPoolSP.dvd1[FC_Use]+MeanforPoolSP.dvd2[FC_Use])/2
	MeanforPool=(MeanforPoolSP.dvd1+MeanforPoolSP.dvd2)/2
	FC_Use2=which(Var_FC_Use>=Mean_FC_Use)
	Var_FC_Use2=Var_FC_Use[FC_Use2]
	Mean_FC_Use2=Mean_FC_Use[FC_Use2]
	Phi=mean((Var_FC_Use2-Mean_FC_Use2)/Mean_FC_Use2^2)
	VarEst=	MeanforPool*(1+MeanforPool*Phi)
	if(Print==T)message(paste("No Replicate - estimate phi",round(Phi,5), "\n"))
	names(VarEst)=names(MeanforPoolSP.dvd1)=
	names(MeanforPoolSP.dvd2)=rownames(DataList.unlist.dvd)
	}

	#DataListSP Here also unlist.. Only two lists
	DataListSP=vector("list",nlevels(Conditions))
	DataListSP.dvd=vector("list",nlevels(Conditions))
	SizeFSP=DataListSP
	MeanSP=DataListSP
	VarSP=DataListSP
	GetPSP=DataListSP
	RSP=DataListSP
	CISP=DataListSP
	tauSP=DataListSP
	NumSampleEachCon=rep(NULL,nlevels(Conditions))

	for (lv in 1:nlevels(Conditions)){
		DataListSP[[lv]]= matrix(DataList.unlist[,Conditions==levels(Conditions)[lv]],nrow=dim(DataList.unlist)[1])
		rownames(DataListSP[[lv]])=rownames(DataList.unlist)
		DataListSP.dvd[[lv]]= matrix(DataList.unlist.dvd[,Conditions==levels(Conditions)[lv]],nrow=dim(DataList.unlist.dvd)[1])
		NumSampleEachCon[lv]=ncol(DataListSP[[lv]])

	if(ncol(DataListSP[[lv]])==1 & !is.null(CI)){
		CISP[[lv]]=matrix(CI[,Conditions==levels(Conditions)[lv]],nrow=dim(DataList.unlist.dvd)[1])
		tauSP[[lv]]=matrix(tau[,Conditions==levels(Conditions)[lv]],nrow=dim(DataList.unlist.dvd)[1])
	}
	# no matter sizeFactors is a vector or a matrix. Matrix should be columns are the normalization factors
	# may input one for each 
	if(length(sizeFactors)==ncol(Data))SizeFSP[[lv]]=sizeFactors[Conditions==levels(Conditions)[lv]]
	if(length(sizeFactors)==length(Data))SizeFSP[[lv]]=sizeFactors[,Conditions==levels(Conditions)[lv]]
	
	
	MeanSP[[lv]]=rowMeans(DataListSP.dvd[[lv]])
	names(MeanSP[[lv]])=rownames(DataListSP[[lv]])
		
	if(length(sizeFactors)==ncol(Data))PrePareVar=sapply(1:ncol( DataListSP[[lv]]),function(i)( DataListSP[[lv]][,i]- SizeFSP[[lv]][i]*MeanSP[[lv]])^2 /SizeFSP[[lv]][i])
	if(length(sizeFactors)==length(Data))PrePareVar=sapply(1:ncol( DataListSP[[lv]]),function(i)( DataListSP[[lv]][,i]- SizeFSP[[lv]][,i]*MeanSP[[lv]])^2 /SizeFSP[[lv]][,i])

	if(ncol(DataListSP[[lv]])==1 & !is.null(CI))
		VarSP[[lv]]=as.vector(((DataListSP[[lv]]/tauSP[[lv]]) * CISP[[lv]]/(CIthre*2))^2)
	if(ncol(DataListSP[[lv]])!=1){
		VarSP[[lv]]=rowSums(PrePareVar)/ncol( DataListSP[[lv]])
		names(MeanSP[[lv]])=rownames(DataList.unlist)
		names(VarSP[[lv]])=rownames(DataList.unlist)
		GetPSP[[lv]]=MeanSP[[lv]]/VarSP[[lv]]
		RSP[[lv]]=MeanSP[[lv]]*GetPSP[[lv]]/(1-GetPSP[[lv]])
	}
}
	
	
	VarList=apply(DataList.unlist.dvd, 1, var)
	if(ncol(Data)==2){
		PoolVar=VarEst
		VarSP[[1]]=VarSP[[2]]=VarEst
		GetPSP[[1]]=MeanSP[[1]]/VarEst
		GetPSP[[2]]=MeanSP[[2]]/VarEst

	}
	if(!ncol(Data)==2){
		CondWithRep=which(NumSampleEachCon>1)
		VarCondWithRep=do.call(cbind,VarSP[CondWithRep])
		PoolVar=rowMeans(VarCondWithRep)
	
	}
	GetP=MeanList/PoolVar
	
    EmpiricalRList=MeanList*GetP/(1-GetP) 
	EmpiricalRList[EmpiricalRList==Inf]	=max(EmpiricalRList[EmpiricalRList!=Inf])
#####################
	if(ncol(Data)!=2){
	Varcbind=do.call(cbind,VarSP)
	VarrowMin=apply(Varcbind,1,min)
	}

	if(ncol(Data)==2){
		Varcbind=VarEst
		VarrowMin=VarEst
		VarSP[[1]]=VarSP[[2]]=VarEst
		names(MeanSP[[1]])=names(VarSP[[1]])
		names(MeanSP[[2]])=names(VarSP[[2]])
	}
	# 
	# 
	GoodData=names(MeanList)[EmpiricalRList>0 &  VarrowMin!=0 & EmpiricalRList!=Inf & !is.na(VarrowMin) & !is.na(EmpiricalRList)]
	NotIn=names(MeanList)[EmpiricalRList<=0 | VarrowMin==0 | EmpiricalRList==Inf |  is.na(VarrowMin) | is.na(EmpiricalRList)]
	#print(paste("ZeroVar",sum(VarrowMin==0), "InfR", length(which(EmpiricalRList==Inf)), "Poi", length(which(EmpiricalRList<0)), ""))
	EmpiricalRList.NotIn=EmpiricalRList[NotIn]
	EmpiricalRList.Good=EmpiricalRList[GoodData]
	EmpiricalRList.Good[EmpiricalRList.Good<1]=1+EmpiricalRList.Good[EmpiricalRList.Good<1]
	if(length(sizeFactors)==ncol(Data)){
		EmpiricalRList.Good.mat= outer(EmpiricalRList.Good, sizeFactors)	
		EmpiricalRList.mat= outer(EmpiricalRList, sizeFactors)
	}
	if(length(sizeFactors)==length(Data)){
	EmpiricalRList.Good.mat=EmpiricalRList.Good* sizeFactors[GoodData,]
	EmpiricalRList.mat=EmpiricalRList* sizeFactors
	}

	# Only Use Data has Good q's
	DataList.In=sapply(1:NoneZeroLength, function(i)DataList[[i]][GoodData[GoodData%in%rownames(DataList[[i]])],],simplify=F)
	DataList.NotIn=sapply(1:NoneZeroLength, function(i)DataList[[i]][NotIn[NotIn%in%rownames(DataList[[i]])],],simplify=F)
	DataListIn.unlist=do.call(rbind, DataList.In)
	DataListNotIn.unlist=do.call(rbind, DataList.NotIn)
	
	DataListSPIn=vector("list",nlevels(Conditions))
	DataListSPNotIn=vector("list",nlevels(Conditions))
	EmpiricalRList.Good.mat.SP=EmpiricalRList.mat.SP=vector("list",nlevels(Conditions))
	for (lv in 1:nlevels(Conditions)){
		DataListSPIn[[lv]]= matrix(DataListIn.unlist[,Conditions==levels(Conditions)[lv]],nrow=dim(DataListIn.unlist)[1])
		if(length(NotIn)>0){	
			DataListSPNotIn[[lv]]= matrix(DataListNotIn.unlist[,Conditions==levels(Conditions)[lv]],nrow=dim(DataListNotIn.unlist)[1])
			rownames(DataListSPNotIn[[lv]])=rownames(DataListNotIn.unlist)
	}
		rownames(DataListSPIn[[lv]])=rownames(DataListIn.unlist)
		EmpiricalRList.Good.mat.SP[[lv]]=matrix(EmpiricalRList.Good.mat[,Conditions==levels(Conditions)[lv]],nrow=dim(EmpiricalRList.Good.mat)[1])
		EmpiricalRList.mat.SP[[lv]]=matrix(EmpiricalRList.mat[,Conditions==levels(Conditions)[lv]],nrow=dim(EmpiricalRList.mat)[1])
	}	

	NumOfEachGroupIn=sapply(1:NoneZeroLength, function(i)max(0,dim(DataList.In[[i]])[1]))
	NumOfEachGroupNotIn=sapply(1:NoneZeroLength, function(i)max(0,dim(DataList.NotIn[[i]])[1]))
	

#################
# For output
#################
RealName.EmpiricalRList=sapply(1:NoneZeroLength,function(i)EmpiricalRList[names(EmpiricalRList)%in%NameList[[i]]], simplify=F)
RealName.MeanList=sapply(1:NoneZeroLength,function(i)MeanList[names(MeanList)%in%NameList[[i]]], simplify=F)
RealName.C1MeanList=sapply(1:NoneZeroLength,function(i)MeanSP[[1]][names(MeanSP[[1]])%in%NameList[[i]]], simplify=F)
RealName.C2MeanList=sapply(1:NoneZeroLength,function(i)MeanSP[[2]][names(MeanSP[[2]])%in%NameList[[i]]], simplify=F)
RealName.C1VarList=sapply(1:NoneZeroLength,function(i)VarSP[[1]][names(VarSP[[1]])%in%NameList[[i]]], simplify=F)
RealName.C2VarList=sapply(1:NoneZeroLength,function(i)VarSP[[2]][names(VarSP[[2]])%in%NameList[[i]]], simplify=F)
RealName.DataList=sapply(1:NoneZeroLength,function(i)DataList[[i]][rownames(DataList[[i]])%in%NameList[[i]],], simplify=F)



RealName.VarList=sapply(1:NoneZeroLength,function(i)VarList[names(VarList)%in%NameList[[i]]], simplify=F)
RealName.PoolVarList=sapply(1:NoneZeroLength,function(i)PoolVar[names(PoolVar)%in%NameList[[i]]], simplify=F)


RealName.QList1=sapply(1:NoneZeroLength,function(i)GetPSP[[1]][names(GetPSP[[1]])%in%NameList[[i]]], simplify=F)
RealName.QList2=sapply(1:NoneZeroLength,function(i)GetPSP[[2]][names(GetPSP[[2]])%in%NameList[[i]]], simplify=F)


for (i in 1:NoneZeroLength){
tmp=NameList[[i]]
names=IsoNamesIn[tmp]

RealName.MeanList[[i]]=RealName.MeanList[[i]][NameList[[i]]]
RealName.VarList[[i]]=RealName.VarList[[i]][NameList[[i]]]
RealName.QList1[[i]]=RealName.QList1[[i]][NameList[[i]]]
RealName.QList2[[i]]=RealName.QList2[[i]][NameList[[i]]]
RealName.EmpiricalRList[[i]]=RealName.EmpiricalRList[[i]][NameList[[i]]]
RealName.C1MeanList[[i]]=RealName.C1MeanList[[i]][NameList[[i]]]
RealName.C2MeanList[[i]]=RealName.C2MeanList[[i]][NameList[[i]]]
RealName.PoolVarList[[i]]=RealName.PoolVarList[[i]][NameList[[i]]]
RealName.C1VarList[[i]]=RealName.C1VarList[[i]][NameList[[i]]]
RealName.C2VarList[[i]]=RealName.C2VarList[[i]][NameList[[i]]]
RealName.DataList[[i]]=RealName.DataList[[i]][NameList[[i]],]

names(RealName.MeanList[[i]])=names
names(RealName.VarList[[i]])=names
if(ncol(DataListSP[[1]])!=1){
	names(RealName.QList1[[i]])=names
	names(RealName.C1VarList[[i]])=names
}
if(ncol(DataListSP[[2]])!=1){
	names(RealName.QList2[[i]])=names
	names(RealName.C2VarList[[i]])=names
}

names(RealName.EmpiricalRList[[i]])=names
names(RealName.C1MeanList[[i]])=names
names(RealName.C2MeanList[[i]])=names
names(RealName.PoolVarList[[i]])=names
rownames(RealName.DataList[[i]])=names


}

#####################
# If Don need EM
#####################	
	if(!is.null(Alpha)&!is.null(Beta)){
	F0Log=f0(Input=DataList.unlist, AlphaIn=Alpha, BetaIn=Beta, 
			 EmpiricalR=EmpiricalRList.mat, NumOfGroups=NumEachGroup, log=T)
	F1Log=f1(Input1=DataListSP[[1]], Input2=DataListSP[[2]], 
			 AlphaIn=Alpha, BetaIn=Beta, EmpiricalRSP1=EmpiricalRList.mat.SP[[1]], 
			 EmpiricalRSP2=EmpiricalRList.mat.SP[[2]], NumOfGroup=NumEachGroup, log=T)
	F0LogMdf=F0Log+600
	F1LogMdf=F1Log+600
	F0Mdf=exp(F0LogMdf)
	F1Mdf=exp(F1LogMdf)
	if(!is.null(PInput)){
		z.list=PInput*F1Mdf/(PInput*F1Mdf+(1-PInput)*F0Mdf)
		PIn=PInput
	}
	if(is.null(PInput)){
		PIn=.5
		PInput=rep(NULL,maxround)
		for(i in 1:maxround){
			z.list=PIn*F1Mdf/(PIn*F1Mdf+(1-PIn)*F0Mdf)
			zNaNName=names(z.list)[is.na(z.list)]
			zGood=which(!is.na(z.list))
			PIn=sum(z.list[zGood])/length(z.list[zGood])
			PInput[i]=PIn
		}
	
	zNaNName=names(z.list)[is.na(z.list)]
	if(length(zNaNName)!=0){
	PNotIn=rep(1-ApproxVal,length(zNaNName))
	MeanList.NotIn=MeanList[zNaNName]
	R.NotIn.raw=MeanList.NotIn*PNotIn/(1-PNotIn)
	if(length(sizeFactors)==ncol(Data))
			R.NotIn=outer(R.NotIn.raw,sizeFactors)
	if(length(sizeFactors)==length(Data))
			R.NotIn=R.NotIn.raw*sizeFactors[zNaNName,]
	R.NotIn1=matrix(R.NotIn[,Conditions==levels(Conditions)[1]],nrow=nrow(R.NotIn))
    R.NotIn2=matrix(R.NotIn[,Conditions==levels(Conditions)[2]],nrow=nrow(R.NotIn))
	NumOfEachGroupNA=sapply(1:NoneZeroLength, function(i)sum(zNaNName%in%rownames(DataList[[i]])))
	F0LogNA=f0(matrix(DataList.unlist[zNaNName,], ncol=ncol(DataList.unlist)), Alpha, Beta, R.NotIn, NumOfEachGroupNA, log=T)
    F1LogNA=f1(matrix(DataListSP[[1]][zNaNName,],ncol=ncol(DataListSP[[1]])), 
               matrix(DataListSP[[2]][zNaNName,],ncol=ncol(DataListSP[[2]])),
			   Alpha, Beta, R.NotIn1,R.NotIn2, NumOfEachGroupNA, log=T)
	F0LogMdfNA=F0LogNA+600
		F1LogMdfNA=F1LogNA+600
		F0MdfNA=exp(F0LogMdfNA)
		F1MdfNA=exp(F1LogMdfNA)
		z.list.NotIn=PIn*F1MdfNA/(PIn*F1MdfNA+(1-PIn)*F0MdfNA)
	z.list[zNaNName]=z.list.NotIn
	F0Log[zNaNName]=F0LogNA
	F1Log[zNaNName]=F1LogNA
	}
	}
	RealName.Z.output=z.list
	RealName.F0=F0Log
	RealName.F1=F1Log
	names(RealName.Z.output)=IsoNamesIn
	names(RealName.F0)=IsoNamesIn
	names(RealName.F1)=IsoNamesIn


	output=list(Alpha=Alpha,Beta=Beta,P=PInput, Z=RealName.Z.output,
			 RList=RealName.EmpiricalRList, MeanList=RealName.MeanList, 
			 VarList=RealName.VarList, QList1=RealName.QList1, QList2=RealName.QList2, 
			 C1Mean=RealName.C1MeanList, C2Mean=RealName.C2MeanList,
			 C1EstVar=RealName.C1VarList, C2EstVar=RealName.C2VarList, 
			 PoolVar=RealName.PoolVarList , DataList=RealName.DataList,
			 PPDE=RealName.Z.output,f0=RealName.F0, f1=RealName.F1)
	return(output)
	}


#####################	
#Initialize SigIn & ...
#####################	
	AlphaIn=0.5
	BetaIn=rep(0.5,NoneZeroLength)
	PIn=0.5


#####################	
# EM
#####################	
	UpdateAlpha=NULL
	UpdateBeta=NULL
	UpdateP=NULL
	UpdatePFromZ=NULL
    Timeperround=NULL 
	for (times in 1:maxround){
    	temptime1=proc.time()
		UpdateOutput=suppressWarnings(LogN(DataListIn.unlist,DataListSPIn, EmpiricalRList.Good.mat ,EmpiricalRList.Good.mat.SP,  NumOfEachGroupIn, AlphaIn, BetaIn, PIn, NoneZeroLength))
    message(paste("iteration", times, "done \n",sep=" "))
		AlphaIn=UpdateOutput$AlphaNew
    	BetaIn=UpdateOutput$BetaNew
    	PIn=UpdateOutput$PNew
		PFromZ=UpdateOutput$PFromZ
    	F0Out=UpdateOutput$F0Out
		F1Out=UpdateOutput$F1Out
		UpdateAlpha=rbind(UpdateAlpha,AlphaIn)
   		UpdateBeta=rbind(UpdateBeta,BetaIn)
    	UpdateP=rbind(UpdateP,PIn)
		UpdatePFromZ=rbind(UpdatePFromZ,PFromZ)
		temptime2=proc.time()
		Timeperround=c(Timeperround,temptime2[3]-temptime1[3])
		message(paste("time" ,round(Timeperround[times],2),"\n",sep=" "))
		Z.output=UpdateOutput$ZNew.list[!is.na(UpdateOutput$ZNew.list)]
   		Z.NA.Names=UpdateOutput$zNaNName
		}
		#Remove this } after testing!!
		 
#    	if (times!=1){  
#        	if((UpdateAlpha[times]-UpdateAlpha[times-1])^2+UpdateBeta[times]-UpdateBeta[times-1])^2+UpdateR[times]-UpdateR[times-1])^2+UpdateP[times]-UpdateP[times-1])^2<=10^(-6)){ 
#           		Result=list(Sig=SigIn, Miu=MiuIn, Tau=TauIn)
#           		break
#        }
#    }
#}

##########Change Names############
## Only z are for Good Ones
GoodData=GoodData[!GoodData%in%Z.NA.Names]
IsoNamesIn.Good=IsoNamesIn[GoodData]
RealName.Z.output=Z.output
RealName.F0=F0Out
RealName.F1=F1Out
names(RealName.Z.output)=IsoNamesIn.Good
names(RealName.F0)=IsoNamesIn.Good
names(RealName.F1)=IsoNamesIn.Good



#########posterior part for other data set here later############
AllNA=unique(c(Z.NA.Names,NotIn))
z.list.NotIn=NULL
AllF0=c(RealName.F0)
AllF1=c(RealName.F1)
AllZ=RealName.Z.output

if (length(AllNA)>0){
	Ng.NA=NgVector[AllNA]
	AllNA.Ngorder=AllNA[order(Ng.NA)]
	NumOfEachGroupNA=rep(0,NoneZeroLength)
	NumOfEachGroupNA.tmp=tapply(Ng.NA,Ng.NA,length)
	names(NumOfEachGroupNA)=c(1:NoneZeroLength)
	NumOfEachGroupNA[names(NumOfEachGroupNA.tmp)]=NumOfEachGroupNA.tmp
	PNotIn=rep(1-ApproxVal,length(AllNA.Ngorder))
	MeanList.NotIn=MeanList[AllNA.Ngorder]
	R.NotIn.raw=MeanList.NotIn*PNotIn/(1-PNotIn) 
	if(length(sizeFactors)==ncol(Data))
	R.NotIn=outer(R.NotIn.raw,sizeFactors)
	if(length(sizeFactors)==length(Data))
	R.NotIn=R.NotIn.raw*sizeFactors[NotIn,]
	R.NotIn1=matrix(R.NotIn[,Conditions==levels(Conditions)[1]],nrow=nrow(R.NotIn))
	R.NotIn2=matrix(R.NotIn[,Conditions==levels(Conditions)[2]],nrow=nrow(R.NotIn))
    
	DataListNotIn.unlistWithZ=matrix(DataList.unlist[AllNA.Ngorder,],nrow=length(AllNA.Ngorder))
	DataListSPNotInWithZ=vector("list",nlevels(Conditions))
	for (lv in 1:nlevels(Conditions)) 
		DataListSPNotInWithZ[[lv]] = matrix(DataListSP[[lv]][AllNA.Ngorder,],nrow=length(AllNA.Ngorder))
		F0Log=f0(DataListNotIn.unlistWithZ,  AlphaIn, BetaIn, R.NotIn, NumOfEachGroupNA, log=T)
    	F1Log=f1(DataListSPNotInWithZ[[1]], DataListSPNotInWithZ[[2]], AlphaIn, BetaIn, R.NotIn1,R.NotIn2, NumOfEachGroupNA, log=T)
		F0LogMdf=F0Log+600
		F1LogMdf=F1Log+600
		F0Mdf=exp(F0LogMdf)
		F1Mdf=exp(F1LogMdf)
		z.list.NotIn=PIn*F1Mdf/(PIn*F1Mdf+(1-PIn)*F0Mdf)
#	names(z.list.NotIn)=IsoNamesIn.Good=IsoNamesIn[which(Names%in%NotIn)]
	names(z.list.NotIn)=IsoNamesIn[AllNA.Ngorder]

	AllZ=c(RealName.Z.output,z.list.NotIn)
	AllZ=AllZ[IsoNamesIn]
	AllZ[is.na(AllZ)]=0
	F0.NotIn=F0Log
	F1.NotIn=F1Log
	names(F0.NotIn)=IsoNamesIn[names(F0Log)]
    names(F1.NotIn)=IsoNamesIn[names(F1Log)]
	AllF0=c(RealName.F0,F0.NotIn)
	AllF1=c(RealName.F1,F1.NotIn)
	AllF0=AllF0[IsoNamesIn]
	AllF1=AllF1[IsoNamesIn]
	AllF0[is.na(AllF0)]=0
	AllF1[is.na(AllF1)]=0
}
PPMatNZ=cbind(1-AllZ,AllZ)
colnames(PPMatNZ)=c("PPEE","PPDE")
rownames(UpdateAlpha)=paste("iter",1:nrow(UpdateAlpha),sep="")
rownames(UpdateBeta)=paste("iter",1:nrow(UpdateBeta),sep="")
rownames(UpdateP)=paste("iter",1:nrow(UpdateP),sep="")
rownames(UpdatePFromZ)=paste("iter",1:nrow(UpdatePFromZ),sep="")
colnames(UpdateBeta)=paste("Ng",1:ncol(UpdateBeta),sep="")

CondOut=levels(Conditions)
names(CondOut)=paste("Condition",c(1:length(CondOut)),sep="")

PPMat=matrix(NA,ncol=2,nrow=nrow(Dataraw))
rownames(PPMat)=rownames(Dataraw)
colnames(PPMat)=c("PPEE","PPDE")
if(is.null(AllZeroNames))PPMat=PPMatNZ
if(!is.null(AllZeroNames))PPMat[names(NotAllZeroNames),]=PPMatNZ[names(NotAllZeroNames),]


#############Result############################
Result=list(Alpha=UpdateAlpha,Beta=UpdateBeta,P=UpdateP,
			PFromZ=UpdatePFromZ, Z=RealName.Z.output,PoissonZ=z.list.NotIn, 
			RList=RealName.EmpiricalRList, MeanList=RealName.MeanList, 
			VarList=RealName.VarList, QList1=RealName.QList1, QList2=RealName.QList2, 
			C1Mean=RealName.C1MeanList, C2Mean=RealName.C2MeanList,C1EstVar=RealName.C1VarList, 
			C2EstVar=RealName.C2VarList, PoolVar=RealName.PoolVarList , 
			DataList=RealName.DataList,PPDE=AllZ,f0=AllF0, f1=AllF1,
			AllZeroIndex=AllZeroNames,PPMat=PPMatNZ, PPMatWith0=PPMat,
			ConditionOrder=CondOut, Conditions=Conditions)
}

