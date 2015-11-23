\name{MedianNorm}
\alias{MedianNorm}
\title{
Median Normalization
}
\description{
'MedianNorm' specifies the median-by-ratio normalization function from Anders et. al., 2010.
}
\usage{
MedianNorm(Data, alternative = FALSE)
}
\arguments{
  \item{Data}{The data matrix with transcripts in rows and lanes in columns.}
	\item{alternative}{if alternative = TRUE, the alternative version of median normalization will be applied. 
	The alternative method is similar to median-by-ratio normalization, but can deal with the cases when all of the
	genes/isoforms have at least one zero counts (in which case the median-by-ratio normalization will fail).
	
	In more details, in median-by-ratio normalization (denote l_1 as libsize for sample
1 as an example, assume total S samples):

\hat{l_1} = median_g [ X_g1 / (X_g1*X_g2*...*X_gS)^{-S} ]     (1)

which estimates l_1 / (l_1 * l_2 * ... * l_S)^{-S}. 
Since we have the constrain that (l_1 * l_2 * ... * l_S) = 1, equation
(1) estimates l_1. Note (1) could also be written as:

\hat{l_1}  = median_g [ (X_g1/X_g1 * X_g1/X_g2 * .... * X_g1/X_gS)^{-S}]

In the alternative method, we estimate l_1/l_1, l_1/l_2, ... l_1/l_S
individually by taking median_g(X_g1/X_g1),  median_g(X_g1/X_g2) ...
Then estimate l_1 = l_1 / (l_1 * l_2 * ... * l_S)^{-S} by taking the
geomean of these estimates:

\hat{l_1} = [ median_g(X_g1/X_g1) * median_g(X_g1/X_g2) *
median_g(X_g1/X_g3) * ... * median_g(X_g1/X_gS) ] ^{-S}

	}
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
