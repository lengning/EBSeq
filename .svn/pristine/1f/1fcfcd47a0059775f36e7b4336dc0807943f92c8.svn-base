QQP <-
function(EBOut, GeneLevel=F){
    if(!"Alpha"%in%names(EBOut))stop("The input doesn't seem like an output from EBTest/EBMultiTest")
		maxround=nrow(EBOut$Alpha)
	  AlphaResult=EBOut$Alpha[maxround,]
    BetaResult=EBOut$Beta[maxround,]
	    
	# Multi
	if(!is.null(EBOut$PPpattern)){
		QList=EBOut$QList
	for(i in 1:length(EBOut$QList)){
		for(j in 1:length(EBOut$QList[[i]])){
			  if(GeneLevel==F)Main=paste("Ig",i,"C",j)
		    if(GeneLevel==T)Main=paste("Gene","C",j)
		tmpSize=length(QList[[i]][[j]][QList[[i]][[j]]<1 & !is.na(QList[[i]][[j]])])
		rdpts=rbeta(tmpSize,AlphaResult,BetaResult[i])
		qqplot(QList[[i]][[j]][QList[[i]][[j]]<1], 
		   rdpts,xlab="estimated q's", ylab="simulated q's from fitted beta",
		   main=Main,
		   xlim=c(0,1),ylim=c(0,1))
	fit=lm(sort(rdpts)~sort(QList[[i]][[j]][QList[[i]][[j]]<1  & !is.na(QList[[i]][[j]])]))
	abline(fit,col="red")
		
	}
	}
	}

	if(is.null(EBOut$PPpattern)){
	for(con in 1:2){
		if(con==1)QList=EBOut$QList1
		if(con==2)QList=EBOut$QList2
		for (i in 1:length(BetaResult)){
			if(GeneLevel==F)Main=paste("Ig",i,"C",con)
	    if(GeneLevel==T)Main=paste("Gene","C",con)
			tmpSize=length(QList[[i]][QList[[i]]<1 & !is.na(QList[[i]])])
			rdpts=rbeta(tmpSize,AlphaResult,BetaResult[i])
			qqplot(QList[[i]][QList[[i]]<1], 
	    rdpts,xlab="estimated q's", ylab="simulated q's from fitted beta",
		  main=Main,
		  xlim=c(0,1),ylim=c(0,1))
	fit=lm(sort(rdpts)~sort(QList[[i]][QList[[i]]<1  & !is.na(QList[[i]])]))
	abline(fit,col="red")
	}	}}
}

