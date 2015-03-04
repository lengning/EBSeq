\name{MedianNorm}
\alias{MedianNorm}
\title{
Median Normalization
}
\description{
'MedianNorm' specifies the median normalization function from Anders et. al., 2010.
}
\usage{
MedianNorm(Data)
}
\arguments{
  \item{Data}{The data matrix with transcripts in rows and lanes in columns.}
}

\value{The function will return a vector contains the normalization factor for each lane.}

\references{
Simon Anders and Wolfgang Huber. Differential expression analysis for sequence count data.
Genome Biology (2010) 11:R106 (open access)
}
\author{
Ning Leng
}


\seealso{
QuantileNorm
}
\examples{
data(GeneMat)
Sizes = MedianNorm(GeneMat)
#EBOut = EBTest(Data = GeneMat,
#	Conditions = as.factor(rep(c("C1","C2"), each=5)),
#	sizeFactors = Sizes, maxround = 5)

}
\keyword{ Normalization }
