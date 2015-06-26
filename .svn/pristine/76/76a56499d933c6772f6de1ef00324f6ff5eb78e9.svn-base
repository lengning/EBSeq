\name{f0}
\alias{f0}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The Prior Predictive Distribution of being EE
}
\description{
'f0' gives the Prior Predictive Distribution of being EE.
}
\usage{
f0(Input, AlphaIn, BetaIn, EmpiricalR, NumOfGroups, log)
}
\arguments{
  \item{Input}{Expression Values.}
  \item{AlphaIn, BetaIn, EmpiricalR}{The parameters estimated from last iteration of EM.}
  \item{NumOfGroups}{How many transcripts within each Ng group.}
  \item{log}{If true, will give the log of the output.}
}
\value{
The function will return the prior predictive distribution values of being EE.
}

\references{
Ning Leng, John A. Dawson, James A. Thomson, Victor Ruotti, Anna I. Rissman, Bart M.G. Smits, Jill D. Haag, Michael N. Gould, Ron M. Stewart, and Christina Kendziorski. EBSeq: An empirical Bayes hierarchical model for inference in RNA-seq experiments. Bioinformatics (2013)
}
\author{
Ning Leng
}

\seealso{
f1
}
\examples{
#
#f0(matrix(rnorm(100,100,1),ncol=10), .5, .6,
#	matrix(rnorm(100,200,1),ncol=10), 100, TRUE)
}
