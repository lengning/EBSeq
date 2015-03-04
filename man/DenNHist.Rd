\name{DenNHist}
\alias{DenNHist}
\title{
Density plot to compare the empirical q's and the simulated q's from the fitted beta distribution.
}
\description{
'DenNHist' gives the density plot that compares the empirical q's and the simulated q's from the fitted beta distribution.
}
\usage{
DenNHist(EBOut, GeneLevel = F)
}
\arguments{
  \item{EBOut}{The output of EBTest or EBMultiTest.}
  \item{GeneLevel}{Indicate whether the results are from data at gene level.}
}

\value{
For data with n1 conditions and n2 uncertainty groups, n1*n2 
plots will be generated. Each plot represents a subset of the data.
The empirical estimation of q's will be represented as blue histograms 
and the density of the fitted beta distribution will be represented as 
the green line.
}

\references{
	Ning Leng, John A. Dawson, James A. Thomson, Victor Ruotti, Anna I. Rissman, Bart M.G. Smits, Jill D. Haag, Michael N. Gould, Ron M. Stewart, and Christina Kendziorski. EBSeq: An empirical Bayes hierarchical model for inference in RNA-seq experiments. Bioinformatics (2013)
}
\author{
Ning Leng
}


\seealso{
beta.mom, QQP, EBTest, EBMultiTest
}
\examples{
data(GeneMat)
GeneMat.small = GeneMat[c(500:1000),]
Sizes = MedianNorm(GeneMat.small)
EBOut = EBTest(Data = GeneMat.small,
	Conditions = as.factor(rep(c("C1","C2"), each=5)),
	sizeFactors = Sizes, maxround = 5)
par(mfrow = c(2,2))
DenNHist(EBOut)
}
\keyword{ beta }
