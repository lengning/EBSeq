\name{LikefunMulti}
\alias{LikefunMulti}
\title{
Likelihood Function of the NB-Beta Model In Multiple Condition Test
}
\description{
'LikefunMulti' specifies the Likelihood Function of the NB-Beta Model In Multiple Condition Test.
}
\usage{
LikefunMulti(ParamPool, InputPool)
}

\arguments{
  \item{ParamPool}{The parameters that will be estimated in EM.}
  \item{InputPool}{The control parameters that will not be estimated in EM.}
}

\value{The function will return the log-likelihood.}

\references{
Ning Leng, John A. Dawson, James A. Thomson, Victor Ruotti, Anna I. Rissman, Bart M.G. Smits, Jill D. Haag, Michael N. Gould, Ron M. Stewart, and Christina Kendziorski. EBSeq: An empirical Bayes hierarchical model for inference in RNA-seq experiments. Bioinformatics (2013)
}

\author{
Ning Leng
}

\examples{
#x1 = c(.6,.7,.3)
#Input = matrix(rnorm(100,100,1),ncol=10)
#RIn = matrix(rnorm(100,200,1),ncol=10)
#InputPool = list(list(Input[,1:5],Input[,6:10]), 
#	Input, cbind(rep(.1, 10), rep(.9,10)), 1, 
#	RIn, list(RIn[,1:5],RIn[,6:10]), 
#	10, rbind(c(1,1),c(1,2)))                    
#LikefunMulti(x1, InputPool)

}


