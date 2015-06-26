\name{GetMultiFC}
\alias{GetMultiFC}
\title{
Calculate the Fold Changes for Multiple Conditions
}
\description{
'GetMultiFC' calculates the Fold Changes for each pair of conditions in a multiple condition study.}
\usage{
GetMultiFC(EBMultiOut, SmallNum = 0.01)
}
\arguments{
  \item{EBMultiOut}{The output of EBMultiTest function.}
\item{SmallNum}{A small number will be added for each transcript in each condition to avoid Inf and NA. 
Default is 0.01.}
}
\details{
Provide the FC (adjusted by the normalization factors) for each pair of comparisons.
A small number will be added for each transcript in each condition to avoid Inf and NA.
Default is set to be 0.01.
}
\value{
\item{FCMat}{The FC of each pair of comparison (adjusted by the normalization factors).}
\item{Log2FCMat}{The log 2 FC of each pair of comparison (adjusted by the normalization factors).}
\item{PostFCMat}{The posterior FC of each pair of comparison.}
\item{Log2PostFCMat}{The log 2 posterior FC of each pair of comparison.}
\item{CondMean}{The mean of each transcript within each condition (adjusted by the normalization factors).}
\item{ConditionOrder}{The condition assignment for C1Mean, C2Mean, etc.}
}


\references{
Ning Leng, John A. Dawson, James A. Thomson, Victor Ruotti, Anna I. Rissman, Bart M.G. Smits, Jill D. Haag, Michael N. Gould, Ron M. Stewart, and Christina Kendziorski. EBSeq: An empirical Bayes hierarchical model for inference in RNA-seq experiments. Bioinformatics (2013)
}

\author{
Ning Leng
}


\seealso{
EBMultiTest, PostFC
}
\examples{
data(MultiGeneMat)
MultiGeneMat.small = MultiGeneMat[201:210,]

Conditions = c("C1","C1","C2","C2","C3","C3")

PosParti = GetPatterns(Conditions)
Parti = PosParti[-3,]

MultiSize = MedianNorm(MultiGeneMat.small)

MultiOut = EBMultiTest(MultiGeneMat.small, 
	NgVector=NULL, Conditions=Conditions, 
	AllParti=Parti, sizeFactors=MultiSize, 
	maxround=5)

MultiFC = GetMultiFC(MultiOut)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Posterior Probability }
