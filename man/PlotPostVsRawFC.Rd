\name{PlotPostVsRawFC}
\alias{PlotPostVsRawFC}
\title{
Plot Posterior FC vs FC
}
\description{
'PlotPostVsRawFC' helps the users visualize the posterior FC vs FC in a two condition study.
}
\usage{
PlotPostVsRawFC(EBOut, FCOut)
}
\arguments{
  \item{EBOut}{
The output of EBMultiTest function.
}
  \item{FCOut}{The output of PostFC function.}

}
\value{
A figure shows fold change vs posterior fold change.
}

\references{
Ning Leng, John A. Dawson, James A. Thomson, Victor Ruotti, Anna I. Rissman, Bart M.G. Smits, Jill D. Haag, Michael N. Gould, Ron M. Stewart, and Christina Kendziorski. EBSeq: An empirical Bayes hierarchical model for inference in RNA-seq experiments. Bioinformatics (2013)
}

\author{
Ning Leng
}

\seealso{
PostFC
}
\examples{
data(GeneMat)
GeneMat.small = GeneMat[c(500:600),]
Sizes = MedianNorm(GeneMat.small)
EBOut = EBTest(Data = GeneMat.small,
	Conditions = as.factor(rep(c("C1","C2"), each=5)),
	sizeFactors = Sizes, maxround = 5)
FC = PostFC(EBOut)
PlotPostVsRawFC(EBOut,FC)

}
\keyword{ Posterior Probability }
