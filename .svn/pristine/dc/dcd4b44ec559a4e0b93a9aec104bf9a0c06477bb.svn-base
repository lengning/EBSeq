\name{EBMultiTest}
\alias{EBMultiTest}
\title{
Using EM algorithm to calculate the posterior probabilities 
of interested patterns in a multiple condition study
}
\description{
'EBMultiTest' is built based on the assumption of NB-Beta Empirical Bayes model. It utilizes the EM algorithm to give the posterior probability of the interested patterns.
}
\usage{
EBMultiTest(Data, NgVector = NULL, Conditions, AllParti = NULL, 
	sizeFactors, maxround, Pool = F, NumBin = 1000, 
	ApproxVal=10^-10, PoolLower=.25, PoolUpper = .75, Print=T,Qtrm=1,QtrmCut=0)
}
\arguments{

  \item{Data}{A data matrix contains expression values for each transcript (gene or isoform level). In which rows should be transcripts and columns should be samples.}
  \item{NgVector}{A vector indicates the uncertainty group assignment of each isoform.
e.g. if we use number of isoforms in the host gene to define the uncertainty groups,
suppose the isoform is in a gene with 2 isoforms, Ng of this isoform should be 2.
The length of this vector should be the same as the number of rows in Data.
If it's gene level data, Ngvector could be left as NULL.}
  \item{Conditions}{A vector indicates the condition in which each sample belongs to. }

  \item{AllParti}{A matrix indicates the interested patterns. Columns shoule be conditions and rows should be patterns. The matrix could be obtained by the GetPatterns function. If AllParti=NULL, all possible patterns will be used.}

  \item{sizeFactors}{The normalization factors. It should be a vector with lane specific numbers (the length of the vector should be the same as the number of samples, with the same order as the columns of Data).}
  \item{maxround}{Number of iterations. The default value is 5.
Users should always check the convergency by looking at the Alpha and
Beta in output. If the hyper-parameter estimations are not converged
in 5 iterations, larger number is suggested.}


\item{Pool}{While working without replicates, user could define the Pool = TRUE in the EBTest function to enable pooling.}
\item{NumBin}{By defining NumBin = 1000, EBSeq will group the genes with similar means together into 1,000 bins.}
\item{PoolLower, PoolUpper}{
With the assumption that only subset of the genes are DE in the data set, we take genes whose FC are in the PoolLower - PoolUpper quantile of the FC's  as the candidate genes (default is 25\%-75\%). 

For each bin, the bin-wise variance estimation is defined as the median of the cross condition variance estimations of the candidate genes within that bin. 

We use the cross condition variance estimations for the candidate genes and the bin-wise variance estimations of the host bin for the non-candidate genes.}

\item{ApproxVal}{The variances of the transcripts with mean < var will be approximated as mean/(1-ApproxVal).}

\item{Print}{Whether print the elapsed-time while running the test.}

\item{Qtrm, QtrmCut}{
Transcripts with Qtrm th quantile < = QtrmCut  will be removed before testing. The default value is Qtrm = 1 and QtrmCut=0.
By default setting, transcripts with all 0's
won't be tested.
}
}


\value{
\item{Alpha}{Fitted parameter alpha of the prior beta distribution. Rows are the values for each iteration.}
\item{Beta}{Fitted parameter beta of the prior beta distribution. Rows are the values for each iteration.}
\item{P, PFromZ}{The bayes estimator of following each pattern of interest. Rows are the values for each iteration.}
\item{Z, PoissonZ}{The Posterior Probability of following each pattern of interest for each transcript. (Maybe not in the same order of input).}
\item{RList}{The fitted values of r for each transcript.}
\item{MeanList}{The mean of each transcript. (across conditions).}
\item{VarList}{The variance of each transcript. (across conditions).}
\item{QList}{The fitted q values of each transcript within each condition.}
\item{SPMean}{The mean of each transcript within each condition (adjusted by the normalization factors).}
\item{SPEstVar}{The estimated variance of each transcript within each condition (adjusted by the normalization factors).}
\item{PoolVar}{The variance of each transcript (The pooled value of within condition EstVar).}
\item{DataList}{A List of data that grouped with Ng and bias.}
\item{PPpattern}{The Posterior Probability of following each pattern (columns) for each transcript (rows). Transcripts with expression 0 for all samples are not shown in this matrix.}
\item{f}{The likelihood of likelihood of prior predictive distribution of being each pattern for each transcript. }
\item{AllParti}{The matrix describe the patterns.}
\item{PPpatternWith0}{The Posterior Probability of following each pattern (columns) for each transcript (rows). Transcripts with expression 0 for all samples are shown in this matrix with PP(any_pattrn)=NA.}
\item{ConditionOrder}{The condition assignment for C1Mean, C2Mean, etc.}
}
\references{
Ning Leng, John A. Dawson, James A. Thomson, Victor Ruotti, Anna I. Rissman, Bart M.G. Smits, Jill D. Haag, Michael N. Gould, Ron M. Stewart, and Christina Kendziorski. EBSeq: An empirical Bayes hierarchical model for inference in RNA-seq experiments. Bioinformatics (2013)
}
\author{
Ning Leng
}

\seealso{
EBTest, GetMultiPP, GetMultiFC
}
\examples{
data(MultiGeneMat)
MultiGeneMat.small = MultiGeneMat[201:210,]
Conditions = c("C1","C1","C2","C2","C3","C3")
PosParti = GetPatterns(Conditions)
Parti = PosParti[-3,]
MultiSize = MedianNorm(MultiGeneMat.small)
MultiOut = EBMultiTest(MultiGeneMat.small, NgVector = NULL,
	Conditions = Conditions, AllParti = Parti, 
	sizeFactors = MultiSize, maxround = 5)
MultiPP = GetMultiPP(MultiOut)
}
\keyword{ DE }
\keyword{ Multiple Condition }% __ONLY ONE__ keyword per line
