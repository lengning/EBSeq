\name{QuantileNorm}
\alias{QuantileNorm}
\title{
Quantile Normalization
}
\description{
'QuantileNorm' gives the quantile normalization.
}
\usage{
QuantileNorm(Data, Quantile)
}
\arguments{

  \item{Data}{
The data matrix with transcripts in rows and lanes in columns.
}
\item{Quantile}{
The quantile the user wishs to use. Should be a number between 0 and 1.	
}
}
\details{
Use a quantile point to normalize the data.
}
\value{
The function will return a vector contains the normalization factor for each lane.
% ...
}
\references{
Bullard, James H., et al. Evaluation of statistical methods for normalization and differential expression in mRNA-Seq experiments. BMC bioinformatics 11.1 (2010): 94.
}
\author{
Ning Leng
}


\seealso{
MedianNorm
}
\examples{
data(GeneMat)
Sizes = QuantileNorm(GeneMat,.75)
#EBOut = EBTest(Data = GeneMat,
#	Conditions = as.factor(rep(c("C1","C2"), each=5)),
#	sizeFactors = Sizes, maxround = 5)

}
\keyword{ Normalization }% __ONLY ONE__ keyword per line
