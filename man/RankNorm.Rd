\name{RankNorm}
\alias{RankNorm}
\title{
Rank Normalization
}
\description{
'RankNorm' gives the rank normalization.
}
\usage{
RankNorm(Data)
}
\arguments{
  \item{Data}{
The data matrix with transcripts in rows and lanes in columns.
}
}

\value{
The function will return a matrix contains the normalization factor for each lane and each transcript.
}

\author{
Ning Leng
}


\seealso{
MedianNorm, QuantileNorm
}

\examples{
data(GeneMat)
Sizes = RankNorm(GeneMat)
# Run EBSeq
# EBres = EBTest(Data = GeneData, NgVector = rep(1,10^4), 
#	Vect5End = rep(1,10^4), Vect3End = rep(1,10^4), 
#	Conditions = as.factor(rep(c(1,2), each=5)), 
#	sizeFactors = Sizes, maxround=5)
}
\keyword{ Normalization }
