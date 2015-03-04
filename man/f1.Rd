\name{f1}
\alias{f1}
\title{
The Prior Predictive Distribution of being DE
}
\description{
'f1' gives the Prior Predictive Distribution of DE.
}
\usage{
f1(Input1, Input2, AlphaIn, BetaIn, EmpiricalRSP1, 
	EmpiricalRSP2, NumOfGroup, log)
}
\arguments{

  \item{Input1}{Expressions from Condition1.}
  \item{Input2}{Expressions from Condition2.}
  \item{AlphaIn, BetaIn, EmpiricalRSP1, EmpiricalRSP2}{The parameters estimated from last iteration of EM.}
  \item{NumOfGroup}{ How many transcripts within each Ng group.}
  \item{log}{If true, will give the log of the output.}
}

\value{
The function will return the prior predictive distribution values of being DE.
}

\references{
Ning Leng, John A. Dawson, James A. Thomson, Victor Ruotti, Anna I. Rissman, Bart M.G. Smits, Jill D. Haag, Michael N. Gould, Ron M. Stewart, and Christina Kendziorski. EBSeq: An empirical Bayes hierarchical model for inference in RNA-seq experiments. Bioinformatics (2013)
}
\author{
Ning Leng
}

\seealso{
f0
}
\examples{
#f1(matrix(rnorm(100,100,1),ncol=10), 
#	matrix(rnorm(100,100,1),ncol=10), .5, .6, 
#	matrix(rnorm(100,200,1),ncol=10), 
#	matrix(rnorm(100,200,1),ncol=10), 100, TRUE)

}
