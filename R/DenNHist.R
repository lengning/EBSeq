DenNHist <-
function(EBOut, GeneLevel=F)
{	
	if(!"Alpha"%in%names(EBOut))stop("The input doesn't seem like an output from EBTest/EBMultiTest")
	maxround=nrow(EBOut$Alpha)
	Alpha=EBOut$Alpha[maxround,]
	Beta=EBOut$Beta[maxround,]
	# Multi
	if(!is.null(EBOut$PPpattern)){
		QList=EBOut$QList
	for(i in 1:length(EBOut$QList)){
		for(j in 1:length(EBOut$QList[[i]])){
		if(GeneLevel==F)Main=paste("Ig",i,"C",j)
		if(GeneLevel==T)Main=paste("Gene","C",j)
		hist(QList[[i]][[j]][QList[[i]][[j]]<.98&QList[[i]][[j]]>0],
		prob=T,col="blue",breaks=100,
		main=Main,
		xlim=c(0,1),xlab=paste("Q alpha=",round(Alpha,2),
								" beta=",round(Beta[i],2),sep=""))
 		tmpSize=length(QList[[i]][[j]][QList[[i]][[j]]<.98])
    tmpseq=seq(0.001,1,length=1000)
	  ll=tmpseq
		lines(ll,dbeta(ll,Alpha,Beta[i]),col="green",lwd=2)
	legend("topright",c("Data","Fitted density"),col=c("blue","green"),lwd=2)
		}
	}
	}

	if(is.null(EBOut$PPpattern)){
	for(con in 1:2){
	if(con==1)QList=EBOut$QList1
	if(con==2)QList=EBOut$QList2
  if(!is.list(QList)) QList=list(QList) 	
	for (i in 1:length(QList)){
	    if(GeneLevel==F)Main=paste("Ig",i,"C",con)
	    if(GeneLevel==T)Main=paste("Gene","C",con)
	hist(QList[[i]][QList[[i]]<.98&QList[[i]]>0],
		 prob=T,col="blue",breaks=100,
		 main=Main,
		 xlim=c(0,1),xlab=paste("Q alpha=",round(Alpha,2),
								" beta=",round(Beta[i],2),sep=""))
	tmpSize=length(QList[[i]][QList[[i]]<.98])
  tmpseq=seq(0.001,1,length=1000)
	ll=tmpseq
	lines(ll,dbeta(ll,Alpha,Beta[i]),col="green",lwd=2)
	legend("topright",c("Data","Fitted density"),col=c("blue","green"),lwd=2)
	}}
	}
	
	}

