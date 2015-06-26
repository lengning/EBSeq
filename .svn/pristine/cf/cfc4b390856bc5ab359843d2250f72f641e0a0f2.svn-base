PlotPostVsRawFC<-function(EBOut,FCOut){
library(gplots)

par(fig=c(0,.8,0,1), new=F)
RainbowColor=rev(redgreen(length(FCOut$PostFC)))
par(oma=c(0,1,1,0),cex=1.3)
plot(FCOut$PostFC,FCOut$RealFC,
log="xy",col=RainbowColor[rank(unlist(EBOut$MeanList)[names(FCOut$PostFC)])], xlab="Posterior FC",
ylab="FC",pch=21)
abline(h=1, v=1, col="gray")
#legend("topleft",col=c("green","black","red"),legend=c("Low Expression","Median Expression","High Expression"),pch=21)

par(fig=c(.7,1,0,1), new=TRUE)
Seq=1:ceiling(length(RainbowColor)/100)
plot(c(0,10), c(1,length(RainbowColor)), type='n', bty='n', xaxt='n', xlab='Rank', ylab='')
for (i in 1:length(Seq)) {
rect(0,(i-1)*100,10,i*100, col=RainbowColor[(i-1)*100], border=NA)
}


}
