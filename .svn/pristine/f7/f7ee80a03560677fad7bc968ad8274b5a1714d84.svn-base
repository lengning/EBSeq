\name{LogN}
\alias{LogN}
\title{
The function to run EM (one round) algorithm for the NB-beta model.
}
\description{
'LogN' specifies the function to run (one round of) the EM algorithm for the NB-beta model.
}
\usage{
LogN(Input, InputSP, EmpiricalR, EmpiricalRSP, NumOfEachGroup, 
	AlphaIn, BetaIn, PIn, NoneZeroLength)
}
\arguments{
  \item{Input, InputSP}{The expressions among all the samples.}
  \item{NumOfEachGroup}{Number of genes in each Ng group.}
  \item{AlphaIn, PIn, BetaIn, EmpiricalR, EmpiricalRSP}{The parameters from the last EM step.}
  \item{NoneZeroLength}{Number of Ng groups.}
}

\references{
Ning Leng, John A. Dawson, James A. Thomson, Victor Ruotti, Anna I. Rissman, Bart M.G. Smits, Jill D. Haag, Michael N. Gould, Ron M. Stewart, and Christina Kendziorski. EBSeq: An empirical Bayes hierarchical model for inference in RNA-seq experiments. Bioinformatics (2013)

}
\author{
Ning Leng
}


\examples{

#Input = matrix(rnorm(100,100,1), ncol=10)
#rownames(Input) = paste("g",1:10)
#RIn = matrix(rnorm(100,200,1), ncol=10)
#res = LogN(Input, list(Input[,1:5], Input[,6:10]),
#	RIn, list(RIn[,1:5], RIn[,6:10]), 
#	10, .6, .7, .3, 1)

}

