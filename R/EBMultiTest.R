EBMultiTest <-
function(Data,NgVector=NULL,Conditions,AllParti=NULL, sizeFactors, maxround,  Pool=F, NumBin=1000, ApproxVal=10^-10,PoolLower=.25, PoolUpper=.75,Print=T,Qtrm=1,QtrmCut=0)
{
 expect_is(sizeFactors, c("numeric","integer"))
 expect_is(maxround, "integer")
 if(!is.factor(Conditions))Conditions=as.factor(Conditions)
 if(is.null(rownames(Data)))stop("Please add gene/isoform names to the data matrix")
 if(!is.matrix(Data))stop("The input Data is not a matrix")
 if(length(Conditions)!=ncol(Data))stop("The number of conditions is not the same as the number of samples! ")
 if(nlevels(Conditions)==2)stop("Only 2 conditions - Please use EBTest() function")
 if(nlevels(Conditions)<2)stop("Less than 2 conditions - Please check your input")
 if(length(sizeFactors)!=length(Data) &  length(sizeFactors)!=ncol(Data))
 stop("The number of library size factors is not the same as the number of samples!")

	tau=CI=CIthre=NULL
	Dataraw=Data

	#Normalized
	DataNorm=GetNormalizedMat(Data, sizeFactors)

	QuantileFor0=apply(DataNorm,1,function(i)quantile(i,Qtrm))
	AllZeroNames=which(QuantileFor0<=QtrmCut)
	NotAllZeroNames=which(QuantileFor0>QtrmCut)
	if(length(AllZeroNames)>0 & Print==T)
					    cat(paste0("Removing transcripts with ",Qtrm*100,
							    " th quantile < = ",QtrmCut," \n",
									length(NotAllZeroNames)," transcripts will be tested \n"))
	if(length(NotAllZeroNames)==0)stop("0 transcript passed")
	Data=Data[NotAllZeroNames,]
	
	if(!is.null(NgVector))NgVector=NgVector[NotAllZeroNames]
	if(is.null(NgVector))NgVector=rep(1,nrow(Data))
	if(length(sizeFactors)!=ncol(Data))sizeFactors=sizeFactors[NotAllZeroNames,]

	#ReNameThem
	IsoNamesIn=rownames(Data)
	Names=paste("I",c(1:dim(Data)[1]),sep="")
	names(IsoNamesIn)=Names
	rownames(Data)=paste("I",c(1:dim(Data)[1]),sep="")
	names(NgVector)=paste("I",c(1:dim(Data)[1]),sep="")
	
	# If PossibleCond==NULL, use all combinations
	NumCond=nlevels(Conditions)
	CondLevels=levels(Conditions)
	#library(blockmodeling)
	if(is.null(AllParti)){
		AllPartiList=sapply(1:NumCond,function(i)nkpartitions(NumCond,i))
		AllParti=do.call(rbind,AllPartiList)
		colnames(AllParti)=CondLevels
	    rownames(AllParti)=paste("Pattern",1:nrow(AllParti),sep="")
	}
	if(length(sizeFactors)==length(Data)){
		rownames(sizeFactors)=rownames(Data)
		colnames(sizeFactors)=Conditions
	}

	
	NoneZeroLength=nlevels(as.factor(NgVector))
	NameList=sapply(1:NoneZeroLength,function(i)names(NgVector)[NgVector==i],simplify=F)
	DataList=sapply(1:NoneZeroLength , function(i) Data[NameList[[i]],],simplify=F)
	names(DataList)=names(NameList)
    
	NumEachGroup=sapply(1:NoneZeroLength , function(i)dim(DataList)[i])
	# Unlist 
	DataList.unlist=do.call(rbind, DataList)

	# Divide by SampleSize factor
	
	if(length(sizeFactors)==ncol(Data))
	DataList.unlist.dvd=t(t( DataList.unlist)/sizeFactors)
	
	if(length(sizeFactors)==length(Data))
	DataList.unlist.dvd=DataList.unlist/sizeFactors
	
	# Pool or Not
	if(Pool==T){
	DataforPoolSP.dvd=MeanforPoolSP.dvd=vector("list",NumCond)
	for(lv in 1:NumCond){
		DataforPoolSP.dvd[[lv]]=matrix(DataList.unlist.dvd[,Conditions==levels(Conditions)[lv]],nrow=dim(DataList.unlist)[1])	
		MeanforPoolSP.dvd[[lv]]=rowMeans(DataforPoolSP.dvd[[lv]])
	}
	MeanforPool.dvd=rowMeans(DataList.unlist.dvd)
	NumInBin=floor(dim(DataList.unlist)[1]/NumBin)
	StartSeq=c(0:(NumBin-1))*NumInBin+1
	EndSeq=c(StartSeq[-1]-1,dim(DataList.unlist)[1])
	MeanforPool.dvd.Sort=sort(MeanforPool.dvd,decreasing=T)	
	MeanforPool.dvd.Order=order(MeanforPool.dvd,decreasing=T)
	PoolGroups=sapply(1:NumBin,function(i)(names(MeanforPool.dvd.Sort)[StartSeq[i]:EndSeq[i]]),simplify=F)
	#FCforPool=MeanforPoolSP.dvd1/MeanforPoolSP.dvd2
	# Use GeoMean of every two-group partition
	Parti2=nkpartitions(NumCond,2)
	FCForPoolList=sapply(1:nrow(Parti2),function(i)rowMeans(do.call(cbind,
							MeanforPoolSP.dvd[Parti2[i,]==1]))/
							rowMeans(do.call(cbind,MeanforPoolSP.dvd[Parti2[i,]==2])),
							simplify=F)
	FCForPoolMat=do.call(cbind,FCForPoolList)
	FCforPool=apply(FCForPoolMat,1,function(i)exp(mean(log(i))))
	names(FCforPool)=names(MeanforPool.dvd)
	FC_Use=names(FCforPool)[which(FCforPool>=quantile(FCforPool[!is.na(FCforPool)],PoolLower) & FCforPool<=quantile(FCforPool[!is.na(FCforPool)],PoolUpper))]
	PoolGroupVar=sapply(1:NumBin,function(i)(mean(apply(matrix(DataList.unlist[PoolGroups[[i]][PoolGroups[[i]]%in%FC_Use],],ncol=ncol(DataList.unlist)),1,var))))	
	PoolGroupVarInList=sapply(1:NumBin,function(i)(rep(PoolGroupVar[i],length(PoolGroups[[i]]))),simplify=F)
	PoolGroupVarVector=unlist(PoolGroupVarInList)
	VarPool=PoolGroupVarVector[MeanforPool.dvd.Order]
	names(VarPool)=names(MeanforPool.dvd)
		}

	DataListSP=vector("list",nlevels(Conditions))
	DataListSP.dvd=vector("list",nlevels(Conditions))
	SizeFSP=DataListSP
	MeanSP=DataListSP
	VarSP=DataListSP
	GetPSP=DataListSP
	RSP=DataListSP
	CISP=DataListSP
	tauSP=DataListSP
	
	NumEachCondLevel=summary(Conditions)
	if(Pool==F & is.null(CI)) CondLevelsUse=CondLevels[NumEachCondLevel>1]
	if(Pool==T | !is.null(CI)) CondLevelsUse=CondLevels
	NumCondUse=length(CondLevelsUse)	

	for (lv in 1:nlevels(Conditions)){
	DataListSP[[lv]]= matrix(DataList.unlist[,Conditions==levels(Conditions)[lv]],nrow=dim(DataList.unlist)[1])
	rownames(DataListSP[[lv]])=rownames(DataList.unlist)
	DataListSP.dvd[[lv]]= matrix(DataList.unlist.dvd[,Conditions==levels(Conditions)[lv]],nrow=dim(DataList.unlist.dvd)[1])
	if(ncol(DataListSP[[lv]])==1 & Pool==F & !is.null(CI)){
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

	if(ncol(DataListSP[[lv]])==1 & Pool==F & !is.null(CI))
		VarSP[[lv]]=as.vector(((DataListSP[[lv]]/tauSP[[lv]]) * CISP[[lv]]/(CIthre*2))^2)
	if( Pool==T){
		VarSP[[lv]]=VarPool	
		}
	if(ncol(DataListSP[[lv]])!=1){
		VarSP[[lv]]=rowSums(PrePareVar)/ncol( DataListSP[[lv]])
		names(VarSP[[lv]])=rownames(DataList.unlist)
		GetPSP[[lv]]=MeanSP[[lv]]/VarSP[[lv]]
	    RSP[[lv]]=MeanSP[[lv]]*GetPSP[[lv]]/(1-GetPSP[[lv]])
	}
	names(MeanSP[[lv]])=rownames(DataList.unlist)
	}

	# Get Empirical R
	# POOL R???
	MeanList=rowMeans(DataList.unlist.dvd)
	VarList=apply(DataList.unlist.dvd, 1, var)

	if(NumCondUse!=0){
		Varcbind=do.call(cbind,VarSP[CondLevels%in%CondLevelsUse])
	  PoolVarSpeedUp_MDFPoi_NoNormVarList=rowMeans(Varcbind)
  	VarrowMin=apply(Varcbind,1,min)

	}	
	if(NumCondUse==0)
	{
	NumFCgp=choose(NumCond,2)
	FC_Use_tmp=vector("list",NumFCgp)
	aa=1
	for(k1 in 1:(NumCond-1)){
		for(k2 in (k1+1):NumCond){
	FCforPool=DataList.unlist.dvd[,k1]/DataList.unlist.dvd[,k2]
	names(FCforPool)=rownames(DataList.unlist.dvd)
	FC_Use_tmp[[aa]]=names(FCforPool)[which(FCforPool>=quantile(FCforPool[!is.na(FCforPool)],.25) & 
								    FCforPool<=quantile(FCforPool[!is.na(FCforPool)],.75))]
	aa=aa+1
	}}
	FC_Use=Reduce(intersect,FC_Use_tmp)
	if(length(FC_Use)==0){
		All_candi=unlist(FC_Use_tmp)
		FC_Use=names(table(All_candi))[1:3]

	}
	Var_FC_Use=apply( DataList.unlist.dvd[FC_Use,],1,var )
	MeanforPool=apply( DataList.unlist.dvd,1,mean )
	Mean_FC_Use=apply( DataList.unlist.dvd[FC_Use,],1,mean )
	FC_Use2=which(Var_FC_Use>=Mean_FC_Use)
	Var_FC_Use2=Var_FC_Use[FC_Use2]
	Mean_FC_Use2=Mean_FC_Use[FC_Use2]
	Phi=mean((Var_FC_Use2-Mean_FC_Use2)/Mean_FC_Use2^2)
	VarEst=	MeanforPool*(1+MeanforPool*Phi)
	if(Print==T)message(paste("No Replicate - estimate phi",round(Phi,5), "\n"))
	Varcbind=VarEst
	PoolVarSpeedUp_MDFPoi_NoNormVarList=VarEst
	VarrowMin=VarEst

	}

	
	
	GetP=MeanList/PoolVarSpeedUp_MDFPoi_NoNormVarList
	
    EmpiricalRList=MeanList*GetP/(1-GetP) 
	# sep
	#Rcb=cbind(RSP[[1]],RSP[[2]])
	#Rbest=apply(Rcb,1,function(i)max(i[!is.na(i) & i!=Inf]))
	EmpiricalRList[EmpiricalRList==Inf]	=max(EmpiricalRList[EmpiricalRList!=Inf])
	# fine
	# 
	GoodData=names(MeanList)[EmpiricalRList>0 &  VarrowMin!=0 & EmpiricalRList!=Inf & !is.na(VarrowMin) & !is.na(EmpiricalRList)]
	NotIn=names(MeanList)[EmpiricalRList<=0 | VarrowMin==0 | EmpiricalRList==Inf |  is.na(VarrowMin) | is.na(EmpiricalRList)]
	#NotIn.BestR=Rbest[NotIn.raw]
	#NotIn.fix=NotIn.BestR[which(NotIn.BestR>0)]
	#EmpiricalRList[names(NotIn.fix)]=NotIn.fix
	#print(paste("ZeroVar",sum(VarrowMin==0), "InfR", length(which(EmpiricalRList==Inf)), "Poi", length(which(EmpiricalRList<0)), ""))
	#GoodData=c(GoodData.raw,names(NotIn.fix))
	#NotIn=NotIn.raw[!NotIn.raw%in%names(NotIn.fix)]
	EmpiricalRList.NotIn=EmpiricalRList[NotIn]
	EmpiricalRList.Good=EmpiricalRList[GoodData]
	EmpiricalRList.Good[EmpiricalRList.Good<1]=1+EmpiricalRList.Good[EmpiricalRList.Good<1]
	if(length(sizeFactors)==ncol(Data))
	EmpiricalRList.Good.mat= outer(EmpiricalRList.Good, sizeFactors)	
	if(length(sizeFactors)==length(Data))
	EmpiricalRList.Good.mat=EmpiricalRList.Good* sizeFactors[GoodData,]


	# Only Use Data has Good q's
	DataList.In=sapply(1:NoneZeroLength, function(i)DataList[[i]][GoodData[GoodData%in%rownames(DataList[[i]])],],simplify=F)
	DataList.NotIn=sapply(1:NoneZeroLength, function(i)DataList[[i]][NotIn[NotIn%in%rownames(DataList[[i]])],],simplify=F)
	DataListIn.unlist=do.call(rbind, DataList.In)
	DataListNotIn.unlist=do.call(rbind, DataList.NotIn)
	
	DataListSPIn=vector("list",nlevels(Conditions))
	DataListSPNotIn=vector("list",nlevels(Conditions))
	EmpiricalRList.Good.mat.SP=vector("list",nlevels(Conditions))
	for (lv in 1:nlevels(Conditions)){
		DataListSPIn[[lv]]= matrix(DataListIn.unlist[,Conditions==levels(Conditions)[lv]],nrow=dim(DataListIn.unlist)[1])
		if(length(NotIn)>0)	DataListSPNotIn[[lv]]= matrix(DataListNotIn.unlist[,Conditions==levels(Conditions)[lv]],nrow=dim(DataListNotIn.unlist)[1])
		rownames(DataListSPIn[[lv]])=rownames(DataListIn.unlist)
		if(length(NotIn)>0)rownames(DataListSPNotIn[[lv]])=rownames(DataListNotIn.unlist)
		EmpiricalRList.Good.mat.SP[[lv]]=matrix(EmpiricalRList.Good.mat[,Conditions==levels(Conditions)[lv]],nrow=dim(EmpiricalRList.Good.mat)[1])
	}	

	NumOfEachGroupIn=sapply(1:NoneZeroLength, function(i)max(0,dim(DataList.In[[i]])[1]))
	NumOfEachGroupNotIn=sapply(1:NoneZeroLength, function(i)max(0,dim(DataList.NotIn[[i]])[1]))

	#Initialize SigIn & ...
	AlphaIn=0.5
	BetaIn=rep(0.5,NoneZeroLength)
	PIn=rep(1/nrow(AllParti),nrow(AllParti))

	####use while to make an infinity round?
	UpdateAlpha=NULL
	UpdateBeta=NULL
	UpdateP=NULL
	UpdatePFromZ=NULL
    Timeperround=NULL 
	for (times in 1:maxround){
    	temptime1=proc.time()
		UpdateOutput=suppressWarnings(LogNMulti(DataListIn.unlist,DataListSPIn, EmpiricalRList.Good.mat ,EmpiricalRList.Good.mat.SP,  
							   NumOfEachGroupIn, AlphaIn, BetaIn, PIn, NoneZeroLength, AllParti,Conditions))
    	message(paste("iteration", times, "done \n",sep=" "))
		AlphaIn=UpdateOutput$AlphaNew
    	BetaIn=UpdateOutput$BetaNew
    	PIn=UpdateOutput$PNew
		PFromZ=UpdateOutput$PFromZ
    	FOut=UpdateOutput$FGood
		UpdateAlpha=rbind(UpdateAlpha,AlphaIn)
   		UpdateBeta=rbind(UpdateBeta,BetaIn)
    	UpdateP=rbind(UpdateP,PIn)
		UpdatePFromZ=rbind(UpdatePFromZ,PFromZ)
		temptime2=proc.time()
		Timeperround=c(Timeperround,temptime2[3]-temptime1[3])
		message(paste("time" ,round(Timeperround[times],2),"\n",sep=" "))
		Z.output=UpdateOutput$ZEachGood
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
## Others are for ALL Data
GoodData=GoodData[!GoodData%in%Z.NA.Names]
IsoNamesIn.Good=as.vector(IsoNamesIn[GoodData])
RealName.Z.output=Z.output
RealName.F=FOut
rownames(RealName.Z.output)=IsoNamesIn.Good
rownames(RealName.F)=IsoNamesIn.Good

RealName.EmpiricalRList=sapply(1:NoneZeroLength,function(i)EmpiricalRList[names(EmpiricalRList)%in%NameList[[i]]], simplify=F)
RealName.MeanList=sapply(1:NoneZeroLength,function(i)MeanList[names(MeanList)%in%NameList[[i]]], simplify=F)
RealName.SPMeanList=sapply(1:NoneZeroLength,function(i)sapply(1:length(MeanSP), function(j)MeanSP[[j]][names(MeanSP[[j]])%in%NameList[[i]]],simplify=F), simplify=F)
RealName.SPVarList=sapply(1:NoneZeroLength,function(i)sapply(1:length(VarSP), function(j)VarSP[[j]][names(VarSP[[j]])%in%NameList[[i]]],simplify=F), simplify=F)
RealName.DataList=sapply(1:NoneZeroLength,function(i)DataList[[i]][rownames(DataList[[i]])%in%NameList[[i]],], simplify=F)

RealName.VarList=sapply(1:NoneZeroLength,function(i)VarList[names(VarList)%in%NameList[[i]]], simplify=F)
RealName.PoolVarList=sapply(1:NoneZeroLength,function(i)PoolVarSpeedUp_MDFPoi_NoNormVarList[names(PoolVarSpeedUp_MDFPoi_NoNormVarList)%in%NameList[[i]]], simplify=F)
RealName.QList=sapply(1:NoneZeroLength,function(i)sapply(1:length(GetPSP), function(j)GetPSP[[j]][names(GetPSP[[j]])%in%NameList[[i]]],simplify=F), simplify=F)


for (i in 1:NoneZeroLength){
tmp=NameList[[i]]
Names=IsoNamesIn[tmp]
RealName.MeanList[[i]]=RealName.MeanList[[i]][NameList[[i]]]
RealName.VarList[[i]]=RealName.VarList[[i]][NameList[[i]]]
	for(j in 1:NumCond){
		RealName.SPMeanList[[i]][[j]]=RealName.SPMeanList[[i]][[j]][NameList[[i]]]
		if(!is.null(RealName.QList[[i]][[j]])){
			RealName.QList[[i]][[j]]=RealName.QList[[i]][[j]][NameList[[i]]]
			RealName.SPVarList[[i]][[j]]=RealName.SPVarList[[i]][[j]][NameList[[i]]]
			names(RealName.QList[[i]][[j]])=Names
			names(RealName.SPVarList[[i]][[j]])=Names
		}
		names(RealName.SPMeanList[[i]][[j]])=Names
	}
RealName.EmpiricalRList[[i]]=RealName.EmpiricalRList[[i]][NameList[[i]]]
RealName.PoolVarList[[i]]=RealName.PoolVarList[[i]][NameList[[i]]]
RealName.DataList[[i]]=RealName.DataList[[i]][NameList[[i]],]

names(RealName.MeanList[[i]])=Names
names(RealName.VarList[[i]])=Names

names(RealName.EmpiricalRList[[i]])=Names
names(RealName.PoolVarList[[i]])=Names
rownames(RealName.DataList[[i]])=Names

}


#########posterior part for other data set here later############
AllNA=unique(c(Z.NA.Names,NotIn))
AllZ=NULL
AllF=NULL
if(length(AllNA)==0){
	AllZ=RealName.Z.output[IsoNamesIn,]
	AllF=RealName.F[IsoNamesIn,]
}
ZEachNA=NULL
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
	R.NotIn=matrix(outer(R.NotIn.raw,sizeFactors),nrow=length(AllNA.Ngorder))
	if(length(sizeFactors)==length(Data))
	R.NotIn=matrix(R.NotIn.raw*sizeFactors[names(R.NotIn.raw),],nrow=length(AllNA.Ngorder))
    
	DataListNotIn.unlistWithZ=matrix(DataList.unlist[AllNA.Ngorder,],
																	 nrow=length(AllNA.Ngorder))
	rownames(DataListNotIn.unlistWithZ)=AllNA.Ngorder
	DataListSPNotInWithZ=vector("list",nlevels(Conditions))
	RListSPNotInWithZ=vector("list",nlevels(Conditions))
	for (lv in 1:nlevels(Conditions)) {
		DataListSPNotInWithZ[[lv]] = matrix(DataListSP[[lv]][AllNA.Ngorder,],nrow=length(AllNA.Ngorder))
		RListSPNotInWithZ[[lv]]=matrix(R.NotIn[,Conditions==levels(Conditions)[lv]],nrow=length(AllNA.Ngorder))
	}

	FListNA=sapply(1:nrow(AllParti),function(i)sapply(1:nlevels(as.factor(AllParti[i,])),
		        function(j)f0(do.call(cbind, DataListSPNotInWithZ[AllParti[i,]==j]),AlphaIn, BetaIn,
                do.call(cbind,RListSPNotInWithZ[AllParti[i,]==j]), NumOfEachGroupNA, log=T)),
				                       simplify=F)
	for(ii in 1:length(FListNA))
		FListNA[[ii]]=matrix(FListNA[[ii]],nrow=length(AllNA.Ngorder))
	FPartiLogNA=matrix(sapply(FListNA,rowSums),nrow=length(AllNA.Ngorder))
	FMatNA=exp(FPartiLogNA+600)
	
	rownames(FMatNA)=rownames(DataListNotIn.unlistWithZ)
	PMatNA=matrix(rep(1,nrow(DataListNotIn.unlistWithZ)),ncol=1)%*%matrix(PIn,nrow=1)
	FmultiPNA=matrix(FMatNA*PMatNA,nrow=length(AllNA.Ngorder))
  DenomNA=rowSums(FmultiPNA)
	ZEachNA=matrix(apply(FmultiPNA,2,function(i)i/DenomNA),nrow=length(AllNA.Ngorder))

	rownames(ZEachNA)=IsoNamesIn[AllNA.Ngorder]

	AllZ=rbind(RealName.Z.output,ZEachNA)
	AllZ=AllZ[IsoNamesIn,]
	
	F.NotIn=FPartiLogNA
	rownames(F.NotIn)=IsoNamesIn[rownames(FMatNA)]
	AllF=rbind(RealName.F,F.NotIn)
	AllF=AllF[IsoNamesIn,]

}
colnames(AllZ)=rownames(AllParti)
colnames(AllF)=rownames(AllParti)
rownames(UpdateAlpha)=paste("iter",1:nrow(UpdateAlpha),sep="")
rownames(UpdateBeta)=paste("iter",1:nrow(UpdateBeta),sep="")
rownames(UpdateP)=paste("iter",1:nrow(UpdateP),sep="")
rownames(UpdatePFromZ)=paste("iter",1:nrow(UpdatePFromZ),sep="")
colnames(UpdateBeta)=paste("Ng",1:ncol(UpdateBeta),sep="")

CondOut=levels(Conditions)
names(CondOut)=paste("Condition",c(1:length(CondOut)),sep="")

AllZWith0=matrix(NA,ncol=ncol(AllZ),nrow=nrow(Dataraw))
rownames(AllZWith0)=rownames(Dataraw)
colnames(AllZWith0)=colnames(AllZ)
if(is.null(AllZeroNames))AllZWith0=AllZ
if(!is.null(AllZeroNames))AllZWith0[names(NotAllZeroNames),]=AllZ[names(NotAllZeroNames),]

#############Result############################
Result=list(Alpha=UpdateAlpha,Beta=UpdateBeta,P=UpdateP,PFromZ=UpdatePFromZ, 
			Z=RealName.Z.output,PoissonZ=ZEachNA, RList=RealName.EmpiricalRList, MeanList=RealName.MeanList, 
			VarList=RealName.VarList, QList=RealName.QList, SPMean=RealName.SPMeanList, SPEstVar=RealName.SPVarList, 
		   	PoolVar=RealName.PoolVarList , DataList=RealName.DataList,PPpattern=AllZ,f=AllF, AllParti=AllParti,
			PPMat=AllZ,PPMatWith0=AllZWith0, ConditionOrder=CondOut)
}

