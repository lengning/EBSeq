\name{PolyFitPlot}
\alias{PolyFitPlot}
\title{
Fit the mean-var relationship using polynomial regression
}
\description{
'PolyFitPlot' fits the mean-var relationship using polynomial regression.
}

\usage{
PolyFitPlot(X, Y, nterms, xname = "Estimated Mean", 
	yname = "Estimated Var", pdfname = "", 
	xlim =  c(-1,5), ylim = c(-1,7), ChangeXY = F, 
	col = "red")
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
The first group of values want to be fitted by the polynomial regression (e.g Mean of the data).
}
  \item{Y}{
The second group of values want to be fitted by the polynomial regression (e.g. variance of the data). The length of Y should be the same as the length of X.
}
  \item{nterms}{
How many polynomial terms want to be used.
}
  \item{xname}{
Name of the x axis.
}
  \item{yname}{
Name of the y axis.
}
  \item{pdfname}{
Name of the plot.
}
  \item{xlim}{
The x limits of the plot. 
}
  \item{ylim}{
The y limits of the plot.

}
  \item{ChangeXY}{
If ChangeXY is setted to be TRUE, X will be treated as the dependent variable and Y will be treated as the independent one. Default is FALSE.
}
  \item{col}{
Color of the fitted line.
}
}
\value{The PolyFitPlot function provides a smooth scatter plot of two variables and their best fitting line of polynomial regression.
}

\references{
Ning Leng, John A. Dawson, James A. Thomson, Victor Ruotti, Anna I. Rissman, Bart M.G. Smits, Jill D. Haag, Michael N. Gould, Ron M. Stewart, and Christina Kendziorski. EBSeq: An empirical Bayes hierarchical model for inference in RNA-seq experiments. Bioinformatics (2013)
}
\author{
Ning Leng
}


\examples{
data(IsoList)
str(IsoList)
IsoMat = IsoList$IsoMat
IsoNames = IsoList$IsoNames
IsosGeneNames = IsoList$IsosGeneNames
IsoSizes = MedianNorm(IsoMat)
NgList = GetNg(IsoNames, IsosGeneNames)

IsoNgTrun = NgList$IsoformNgTrun
#IsoEBOut = EBTest(Data = IsoMat.small, 
#	NgVector = IsoNgTrun,
#	Conditions = as.factor(rep(c("C1","C2"), each=5)),
#	sizeFactors = IsoSizes, maxround = 5)

#par(mfrow=c(2,2))
#PolyFitValue = vector("list",3)

#for(i in 1:3)
#	PolyFitValue[[i]] = PolyFitPlot(IsoEBOut$C1Mean[[i]],
#		IsoEBOut$C1EstVar[[i]], 5)

#PolyAll = PolyFitPlot(unlist(IsoEBOut$C1Mean), 
#	unlist(IsoEBOut$C1EstVar), 5)

#lines(log10(IsoEBOut$C1Mean[[1]][PolyFitValue[[1]]$sort]),
#	PolyFitValue[[1]]$fit[PolyFitValue[[1]]$sort], 
#	col="yellow", lwd=2)
#lines(log10(IsoEBOut$C1Mean[[2]][PolyFitValue[[2]]$sort]),
#	PolyFitValue[[2]]$fit[PolyFitValue[[2]]$sort], 
#	col="pink", lwd=2)
#lines(log10(IsoEBOut$C1Mean[[3]][PolyFitValue[[3]]$sort]),
#	PolyFitValue[[3]]$fit[PolyFitValue[[3]]$sort], 
#	col="green", lwd=2)

#legend("topleft",c("All Isoforms","Ng = 1","Ng = 2","Ng = 3"),
#	col = c("red","yellow","pink","green"), 
#	lty=1, lwd=3, box.lwd=2)

}
