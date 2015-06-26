\name{GetNormalizedMat}
\alias{GetNormalizedMat}
\title{
Calculate normalized expression matrix
}
\description{
'GetNormalizedMat' calculates the normalized expression matrix. 
(Note: this matrix is only used for visualization etc. EBTes and EBMultiTest request *un-adjusted* expressions and normalization factors.)
}
\usage{
GetNormalizedMat(Data, Sizes)
}
\arguments{

  \item{Data}{The data matrix with transcripts in rows and lanes in columns.}
\item{Sizes}{A vector contains the normalization factor for each lane.}
}
\value{The function will return a normalized matrix.}
\references{
Ning Leng, John A. Dawson, James A. Thomson, Victor Ruotti, Anna I. Rissman, Bart M.G. Smits, Jill D. Haag, Michael N. Gould, Ron M. Stewart, and Christina Kendziorski. EBSeq: An empirical Bayes hierarchical model for inference in RNA-seq experiments. Bioinformatics (2013)
}
\author{
Ning Leng
}

\examples{

data(GeneMat)
str(GeneMat)
Sizes = MedianNorm(GeneMat)
NormData = GetNormalizedMat(GeneMat, Sizes)
}
\keyword{ Normalization }% __ONLY ONE__ keyword per line
