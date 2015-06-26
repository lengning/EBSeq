\name{EBTest}
\alias{EBTest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Using EM algorithm to calculate the posterior probabilities of being DE
}
\description{
Base on the assumption of NB-Beta Empirical Bayes model, the EM algorithm is used to get the posterior probability of being DE.
}
\usage{
EBTest(Data, NgVector = NULL, Conditions, sizeFactors, maxround, 
	Pool = F, NumBin = 1000, ApproxVal = 10^-10, Alpha = NULL, 
	Beta = NULL, PInput = NULL, RInput = NULL, 
	PoolLower = .25, PoolUpper = .75, Print = T, Qtrm = .75,QtrmCut=10)
}
\arguments{

  \item{Data}{A data matrix contains expression values for each transcript (gene or isoform level). In which rows should be transcripts and columns should be samples.}
  \item{NgVector}{A vector indicates the uncertainty group assignment of each isoform. 
e.g. if we use number of isoforms in the host gene to define the uncertainty groups, suppose the isoform is in a gene with 2 isoforms, Ng of this isoform should be 2. The length of this vector should be the same as the number of rows in Data. If it's gene level data, Ngvector could be left as NULL.}
  \item{Conditions}{A factor indicates the condition which each sample belongs to. }
  \item{sizeFactors}{The normalization factors. It should be a vector with lane specific numbers (the length of the vector should be the same as the number of samples, with the same order as the columns of Data).}
  \item{maxround}{Number of iterations. The default value is 5. Users should always check the convergency by looking at the Alpha and Beta in output. If the hyper-parameter estimations are not converged in 5 iterations, larger number is suggested.}
\item{Pool}{While working without replicates, user could define the Pool = TRUE in the EBTest function to enable pooling.}
\item{NumBin}{By defining NumBin = 1000, EBSeq will group the genes with similar means together into 1,000 bins.}
\item{PoolLower, PoolUpper}{
With the assumption that only subset of the genes are DE in the data set, we take genes whose FC are in the PoolLower - PoolUpper quantile of the FC's  as the candidate genes (default is 25\%-75\%). 

For each bin, the bin-wise variance estimation is defined as the median of the cross condition variance estimations of the candidate genes within that bin. 

We use the cross condition variance estimations for the candidate genes and the bin-wise variance estimations of the host bin for the non-candidate genes.
}

\item{ApproxVal}{The variances of the transcripts with mean < var will be approximated as mean/(1-ApproxVal). }

\item{Alpha, Beta, PInput, RInput}{If the parameters are known and the user doesn't want to estimate them from the data, user could specify them here.}
\item{Print}{Whether print the elapsed-time while running the test.}
\item{Qtrm, QtrmCut}{
Transcripts with Qtrm th quantile < = QtrmCut  will be removed before testing. The default value is Qtrm = 0.75 and QtrmCut=10.
By default setting, transcripts that have >75\% of the samples with expression less than 10
won't be tested.
}
}

\details{For each transcript gi within condition, the model assumes:
X_{gis}|mu_{gi} ~ NB (r_{gi0} * l_s, q_{gi})
q_gi|alpha, beta^N_g ~ Beta (alpha, beta^N_g)
In which the l_s is the sizeFactors of samples.

The function will test "H0: q_{gi}^{C1} = q_{gi}^{C2}" and "H1: q_{gi}^{C1} != q_{gi}^{C2}."
}
\value{
\item{Alpha}{Fitted parameter alpha of the prior beta distribution. Rows are the values for each iteration.}
\item{Beta}{Fitted parameter beta of the prior beta distribution. Rows are the values for each iteration.}
\item{P, PFromZ}{The bayes estimator of being DE. Rows are the values for each iteration.}
\item{Z, PoissonZ}{The Posterior Probability of being DE for each transcript(Maybe not in the same order of input). }
\item{RList}{The fitted values of r for each transcript.}
\item{MeanList}{The mean of each transcript (across conditions).}
\item{VarList}{The variance of each transcript (across conditions).}
\item{QListi1}{The fitted q values of each transcript within condition 1.}
\item{QListi2}{The fitted q values of each transcript within condition 2.}
\item{C1Mean}{The mean of each transcript within Condition 1 (adjusted by normalization factors).}
\item{C2Mean}{The mean of each transcript within Condition 2 (adjusted by normalization factors).}
\item{C1EstVar}{The estimated variance of each transcript within Condition 1 (adjusted by normalization factors).}
\item{C2EstVar}{The estimated variance of each transcript within Condition 2 (adjusted by normalization factors).}
\item{PoolVar}{The variance of each transcript (The pooled value of within condition EstVar).}
\item{DataList}{A List of data that grouped with Ng.}
\item{PPDE}{The Posterior Probability of being DE for each transcript (The same order of input).}
\item{f0,f1}{The likelihood of the prior predictive distribution of being EE or DE (in log scale).}
\item{AllZeroIndex}{The transcript with expression 0 for all samples (which are not tested).}
\item{PPMat}{A matrix contains posterior probabilities of being EE (the first column) or DE (the second column). 
Rows are transcripts. 
Transcripts with expression 0 for all samples are not shown in this matrix.}
\item{PPMatWith0}{A matrix contains posterior probabilities of being EE (the first column) or DE (the second column). 
Rows are transcripts. 
Transcripts with expression 0 for all samples are shown as PP(EE) = PP(DE) = NA in this matrix. 
The transcript order is exactly the same as the order of the input data.}
\item{ConditionOrder}{The condition assignment for C1Mean, C2Mean, etc.}
\item{Conditions}{The input conditions.}
}
\references{
Ning Leng, John A. Dawson, James A. Thomson, Victor Ruotti, Anna I. Rissman, Bart M.G. Smits, Jill D. Haag, Michael N. Gould, Ron M. Stewart, and Christina Kendziorski. EBSeq: An empirical Bayes hierarchical model for inference in RNA-seq experiments. Bioinformatics (2013)
}
\author{
Ning Leng
}


\seealso{
EBMultiTest, PostFC, GetPPMat
}
\examples{
data(GeneMat)
str(GeneMat)
GeneMat.small = GeneMat[c(1:10,511:550),]
Sizes = MedianNorm(GeneMat.small)
EBOut = EBTest(Data = GeneMat.small,
	Conditions = as.factor(rep(c("C1","C2"), each = 5)),
	sizeFactors = Sizes, maxround = 5)
PP = GetPPMat(EBOut)
}
\keyword{ DE }
\keyword{ Two condition }% __ONLY ONE__ keyword per line
