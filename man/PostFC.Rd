\name{PostFC}
\alias{PostFC}
\title{
Calculate the posterior fold change for each transcript across conditions
}
\description{
'PostFC' calculates the posterior fold change for each transcript across conditions.
}
\usage{
PostFC(EBoutput, SmallNum = 0.01)
}
\arguments{
  \item{EBoutput}{
The ourput from function EBTest.
}

\item{SmallNum}{A small number will be added for each transcript in each condition to avoid Inf and NA. Default is 0.01.}

}
\value{
Provide both FC and posterior FC across two conditions.
FC is calculated as (MeanC1+SmallNum)/(MeanC2+SmallNum).
And Posterior FC is calculated as:

# Post alpha P_a_C1 = alpha + r_C1 * n_C1

# Post beta P_b_C1 = beta + Mean_C1 * n_C1

# P_q_C1 = P_a_C1 / (P_a_C1 + P_b_C1)

# Post FC = ((1-P_q_C1)/P_q_c1) / ( (1-P_q_c2)/P_q_c2)

\item{PostFC}{The posterior FC across two conditions.}
\item{RealFC}{The FC across two conditions (adjusted by the normalization factors).}
\item{Direction}{The diretion of FC calculation.}

}
\references{
Ning Leng, John A. Dawson, James A. Thomson, Victor Ruotti, Anna I. Rissman, Bart M.G. Smits, Jill D. Haag, Michael N. Gould, Ron M. Stewart, and Christina Kendziorski. EBSeq: An empirical Bayes hierarchical model for inference in RNA-seq experiments. Bioinformatics (2013)
}
\author{
Ning Leng
}


\seealso{
EBTest, GetMultiFC
}
\examples{
data(GeneMat)
GeneMat.small = GeneMat[c(500:550),]
Sizes = MedianNorm(GeneMat.small)
EBOut = EBTest(Data = GeneMat.small,
	Conditions = as.factor(rep(c("C1","C2"), each=5)),
	sizeFactors = Sizes, maxround = 5)
FC=PostFC(EBOut)
}
\keyword{ Fold Change }
