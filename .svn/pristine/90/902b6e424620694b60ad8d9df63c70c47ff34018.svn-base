\name{crit_fun}
\alias{crit_fun}
\title{
Calculate the soft threshold for a target FDR
}
\description{
'crit_fun' calculates the soft threshold for a target FDR.
}
\usage{
crit_fun(PPEE, thre)
}
\arguments{
  \item{PPEE}{The posterior probabilities of being EE.}
  \item{thre}{The target FDR.}

}
\details{
Regarding a target FDR alpha, both hard threshold and soft threshold could be used. 
If the hard threshold is preferred, user could simply take the transcripts with 
PP(DE) greater than (1-alpha). Using the hard threshold, any DE transcript in the 
list is with FDR <= alpha. 

If the soft threshold is preferred, user could take the transcripts with PP(DE) 
greater than crit_fun(PPEE, alpha). Using the soft threshold, the list of 
DE transcripts is with average FDR alpha.
}
\value{
The adjusted FDR threshold of target FDR.
}
\references{
Ning Leng, John A. Dawson, James A. Thomson, Victor Ruotti, Anna I. Rissman, Bart M.G. Smits, Jill D. Haag, Michael N. Gould, Ron M. Stewart, and Christina Kendziorski. EBSeq: An empirical Bayes hierarchical model for inference in RNA-seq experiments. Bioinformatics (2013)
}
\author{
Ning Leng
}

\examples{
data(GeneMat)
GeneMat.small = GeneMat[c(1:10, 500:600),]
Sizes = MedianNorm(GeneMat.small)
EBOut = EBTest(Data = GeneMat.small,
	Conditions = as.factor(rep(c("C1","C2"), each=5)),
	sizeFactors = Sizes, maxround = 5)
PP = GetPPMat(EBOut)
DEfound = rownames(PP)[which(PP[,"PPDE"] >= 0.95)]
str(DEfound)

SoftThre = crit_fun(PP[,"PPEE"], 0.05)
DEfound_soft = rownames(PP)[which(PP[,"PPDE"] >= SoftThre)]
}
\keyword{ FDR }
