\name{GetPP}
\alias{GetPP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate the Posterior Probability of each transcript.
}
\description{
'GetPP' generates the Posterior Probability of being DE of each transcript based on the EBTest output.
}
\usage{
GetPP(EBout)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{EBout}{The output of EBTest function.}
}

\value{The poster probabilities of being DE.
}
\references{
Ning Leng, John A. Dawson, James A. Thomson, Victor Ruotti, Anna I. Rissman, Bart M.G. Smits, Jill D. Haag, Michael N. Gould, Ron M. Stewart, and Christina Kendziorski. EBSeq: An empirical Bayes hierarchical model for inference in RNA-seq experiments. Bioinformatics (2013)
}
\author{
Ning Leng
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{GetPPMat
}
\examples{
data(GeneMat)
GeneMat.small = GeneMat[c(1:10,500:550),]
Sizes = MedianNorm(GeneMat.small)
EBOut = EBTest(Data = GeneMat.small,
	Conditions = as.factor(rep(c("C1","C2"), each=5)),
	sizeFactors = Sizes, maxround = 5)
PPDE = GetPP(EBOut)
str(PPDE)
head(PPDE)
}
\keyword{ Posterior Probability }
