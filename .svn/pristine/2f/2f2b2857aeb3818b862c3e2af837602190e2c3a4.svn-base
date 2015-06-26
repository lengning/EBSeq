GetMultiPP <- function(EBout){
	if(!"PPpattern"%in%names(EBout))stop("The input doesn't seem like an output from EBMultiTest")

	PP=EBout$PPpattern	
	UnderFlow=which(is.na(rowSums(PP)))
	if(length(UnderFlow)!=0)Good=c(1:nrow(PP))[-UnderFlow]
	else Good=c(1:nrow(PP))
	MAP=rep(NA,nrow(PP))
	names(MAP)=rownames(PP)
	MAP[Good]=colnames(PP)[apply(PP[Good,],1,which.max)]
	MAP[UnderFlow]="NoTest"
	AllParti=EBout$AllParti
	out=list(PP=PP, MAP=MAP,Patterns=AllParti)
}
