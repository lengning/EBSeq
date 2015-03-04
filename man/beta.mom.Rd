\name{beta.mom}
\alias{beta.mom}

\title{
Fit the beta distribution by method of moments
}

\description{
'beta.mom' fits the beta distribution by method of moments.
}

\usage{
beta.mom(qs.in)
}

\arguments{
  \item{qs.in}{A vector contains the numbers that are assumed to follow a beta distribution.}
}

\value{
  \item{alpha.hat}{Returns the estimation of alpha.}
  \item{beta.hat}{Returns the estimation of beta.}
}
\references{
Ning Leng, John A. Dawson, James A. Thomson, Victor Ruotti, Anna I. Rissman, Bart M.G. Smits, Jill D. Haag, Michael N. Gould, Ron M. Stewart, and Christina Kendziorski. EBSeq: An empirical Bayes hierarchical model for inference in RNA-seq experiments. Bioinformatics (2013)
}
\author{
Ning Leng
}


\seealso{
DenNHist, DenNHistTable
}
\examples{
#tmp = rbeta(5, 5, 100)
#param = beta.mom(tmp)
}
\keyword{ beta }
