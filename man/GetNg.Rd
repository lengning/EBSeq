\name{GetNg}
\alias{GetNg}
\title{
Ng Vector
}
\description{
'GetNg' generates the Ng vector for the isoform level data. 
(While using the number of isoform in the host gene to define the uncertainty groups.)
}
\usage{
GetNg(IsoformName, GeneName, TrunThre = 3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{IsoformName}{A vector contains the isoform names.}
  \item{GeneName}{The gene names of the isoforms in IsoformNames (Should be in the same order).}
\item{TrunThre}{The number of uncertainty groups the user wish to define. The default is 3.}
}
\value{
  \item{GeneNg}{The number of isoforms that are contained in each gene.
  }
  \item{GeneNgTrun}{The truncated Ng of each gene. (The genes contain more than 3 isoforms are with Ng 3.)
  }
   \item{IsoformNg}{The Ng of each isoform.}
    \item{IsoformNgTrun}{The truncated Ng of each isoform (could be used to define the uncertainty group assignment).}
}

\references{
Ning Leng, John A. Dawson, James A. Thomson, Victor Ruotti, Anna I. Rissman, Bart M.G. Smits, Jill D. Haag, Michael N. Gould, Ron M. Stewart, and Christina Kendziorski. EBSeq: An empirical Bayes hierarchical model for inference in RNA-seq experiments.  Bioinformatics (2013)
}
\author{
Ning Leng
}

\examples{

data(IsoList)

IsoMat = IsoList$IsoMat
IsoNames = IsoList$IsoNames
IsosGeneNames = IsoList$IsosGeneNames
IsoSizes = MedianNorm(IsoMat)
NgList = GetNg(IsoNames, IsosGeneNames)

#IsoNgTrun = NgList$IsoformNgTrun
#IsoEBOut = EBTest(Data = IsoMat, NgVector = IsoNgTrun,
#	Conditions = as.factor(rep(c("C1","C2"), each=5)),
#	sizeFactors = IsoSizes, maxround = 5)

}

\keyword{ Ng }
