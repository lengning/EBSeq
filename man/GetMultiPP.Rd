\name{GetMultiPP}
\alias{GetMultiPP}
\title{
Posterior Probability of Each Transcript
}
\description{
'GetMultiPP' generates the Posterior Probability of being each pattern of each transcript based on the EBMultiTest output.
}
\usage{
GetMultiPP(EBout)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{EBout}{The output of EBMultiTest function.}

}
\value{
\item{PP}{The poster probabilities of being each pattern.}
\item{MAP}{Gives the most likely pattern.}
\item{Patterns}{The Patterns.}
}
\references{
Ning Leng, John A. Dawson, James A. Thomson, Victor Ruotti, Anna I. Rissman, Bart M.G. Smits, Jill D. Haag, Michael N. Gould, Ron M. Stewart, and Christina Kendziorski. EBSeq: An empirical Bayes hierarchical model for inference in RNA-seq experiments. Bioinformatics (2013)
}
\author{
Ning Leng
}


\seealso{GetPPMat}
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
MultiPP = GetMultiPP(MultiOut)
}
\keyword{ Posterior Probability }
