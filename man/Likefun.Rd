\name{Likefun}
\alias{Likefun}
\title{
Likelihood Function of the NB-Beta Model
}
\description{
'Likefun' specifies the Likelihood Function of the NB-Beta Model.
}
\usage{
Likefun(ParamPool, InputPool)
}
\arguments{
  \item{ParamPool}{The parameters that will be estimated in EM.}
  \item{InputPool}{The control parameters that will not be estimated in EM.}
}

\value{The function will return the log-likelihood.
}

\references{
Ning Leng, John A. Dawson, James A. Thomson, Victor Ruotti, Anna I. Rissman, Bart M.G. Smits, Jill D. Haag, Michael N. Gould, Ron M. Stewart, and Christina Kendziorski. EBSeq: An empirical Bayes hierarchical model for inference in RNA-seq experiments. Bioinformatics (2013)
}

\author{
Ning Leng
}


\examples{
#x1 = c(.6,.7,.3)
#Input = matrix(rnorm(100,100,1), ncol=10)
#RIn = matrix(rnorm(100,200,1), ncol=10)
#InputPool = list(Input[,1:5], Input[,6:10], Input,
#	rep(.1,100), 1, RIn, RIn[,1:5], RIn[,6:10], 100) 
#Likefun(x1, InputPool)
}
