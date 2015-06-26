GetNg<- function(IsoformName, GeneName, TrunThre=3){
	if(length(IsoformName)!=length(GeneName))stop("The length of IsoformName is not the same as the length of GeneName")
	GeneNg = tapply(IsoformName, GeneName, length)
	if(max(GeneNg)<TrunThre)stop("The max Ng is less than the TrunThre")
	IsoformNg = GeneNg[GeneName]
	names(IsoformNg) = IsoformName
	GeneNgTrun=GeneNg
	GeneNgTrun[GeneNgTrun>TrunThre]=TrunThre
	IsoformNgTrun=IsoformNg
	IsoformNgTrun[IsoformNgTrun>TrunThre]=TrunThre
	out=list( GeneNg=GeneNg, GeneNgTrun=GeneNgTrun, IsoformNg=IsoformNg, IsoformNgTrun=IsoformNgTrun)
	}
