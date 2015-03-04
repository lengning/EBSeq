GetDEResults<-function(EBPrelim, FDR=0.05, Method="robust", 
											 FDRMethod="hard", Threshold_FC=0.7, 
											 Threshold_FCRatio=0.3, SmallNum=0.01)
{
   if(!"PPDE"%in%names(EBPrelim))stop("The input doesn't seem like an output from EBTest")

	#################
	Conditions = EBPrelim$Conditions
	Levels = levels(as.factor(Conditions))
	PPcut=FDR
	# normalized data
	GeneMat=EBPrelim$DataNorm


  ###Get DEfound by FDRMethod type
  PP=GetPPMat(EBPrelim)
  if(FDRMethod=="hard")
  {DEfound=rownames(PP)[which(PP[,"PPDE"]>=(1-PPcut))]}
  else{SoftThre=crit_fun(PP[,"PPEE"],PPcut)
  DEfound=rownames(PP)[which(PP[,"PPDE"]>=SoftThre)]}
  
	# classic
	if(Method=="classic"){
	Gene_status=rep("EE",dim(GeneMat)[1])
  names(Gene_status)=rownames(GeneMat)
  Gene_status[DEfound]="DE"
  NoTest_genes=rownames(GeneMat)[!(rownames(GeneMat)%in%rownames(PP))]
  Gene_status[NoTest_genes]="Filtered: Low Expression"	
	  
	PPMatWith0=EBPrelim$PPMatWith0
	PPMatWith0[NoTest_genes,]=c(NA,NA)
	
	return(list(DEfound=DEfound,PPMat=PPMatWith0,Status=Gene_status))
	}
	else{
	###Post_Foldchange
  PostFoldChange=PostFC(EBPrelim)
  PPFC=PostFoldChange$PostFC
  
	OldPPFC=PPFC[DEfound]
	OldPPFC[which(OldPPFC>1)]=1/OldPPFC[which(OldPPFC>1)]

	FilterFC=names(OldPPFC)[which(OldPPFC>Threshold_FC)]

  ###New Fold Change
  NewFC1=apply(matrix(GeneMat[DEfound,which(Conditions==Levels[[1]])]+SmallNum,
											nrow=length(DEfound)),1,median)
  NewFC2=apply(matrix(GeneMat[DEfound,which(Conditions==Levels[[2]])]+SmallNum,
											nrow=length(DEfound)),1,median)
	NewFC=NewFC1/NewFC2
	NewFC[which(NewFC>1)]=1/NewFC[which(NewFC>1)]

  ###FC Ratio
  FCRatio=NewFC/OldPPFC
	FCRatio[which(OldPPFC<NewFC)]=1/FCRatio[which(OldPPFC<NewFC)]
  
	FilterFCR=names(FCRatio)[which(FCRatio<Threshold_FCRatio)]

	###Results format Setting
  Gene_status=rep("EE",dim(GeneMat)[1])
  names(Gene_status)=rownames(GeneMat)
  Gene_status[DEfound]="DE"
  NoTest_genes=rownames(GeneMat)[!(rownames(GeneMat)%in%rownames(PP))]
  Gene_status[NoTest_genes]="Filtered: Low Expression"
  ###Filtering
  Filtered_DEfound=setdiff(DEfound, union(FilterFC, FilterFCR))

	PPMatWith0=EBPrelim$PPMatWith0
	NAGenes=union(NoTest_genes,union(FilterFC, FilterFCR))
	PPMatWith0[NAGenes,]=c(NA,NA)
  
	Gene_status[FilterFC]="Filtered: Fold Change"
	Gene_status[FilterFCR]="Filtered: Fold Change Ratio"

  return(list(DEfound=Filtered_DEfound,PPMat=PPMatWith0,Status=Gene_status))
}}
