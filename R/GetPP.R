GetPP <- function(EBout){
	if(!"PPDE"%in%names(EBout))stop("The input doesn't seem like an output from EBTest")

	PP=EBout$PPDE	
}
