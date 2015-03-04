\name{GetPPMat}
\alias{GetPPMat}
\title{
Posterior Probability of Transcripts
}
\description{
'GetPPMat' generates the Posterior Probability of being each pattern of each transcript based on the EBTest output.
}
\usage{
GetPPMat(EBout)
}
\arguments{
  \item{EBout}{The output of EBTest function.}

}
\value{The poster probabilities of being EE (first column) and DE (second column).
}
\references{
Ning Leng, John A. Dawson, James A. Thomson, Victor Ruotti, Anna I. Rissman, Bart M.G. Smits, Jill D. Haag, Michael N. Gould, Ron M. Stewart, and Christina Kendziorski. EBSeq: An empirical Bayes hierarchical model for inference in RNA-seq experiments. Bioinformatics (2013)
}
\author{
Ning Leng
}

\examples{

data(GeneMat)
GeneMat.small = GeneMat[c(500:550),]
Sizes = MedianNorm(GeneMat.small)
EBOut = EBTest(Data = GeneMat.small,
	Conditions = as.factor(rep(c("C1","C2"), each=5)),
	sizeFactors = Sizes, maxround = 5)
PP = GetPPMat(EBOut)
str(PP)
head(PP)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Posterior Probability }
