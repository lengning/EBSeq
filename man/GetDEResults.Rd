\name{GetDEResults}
\alias{GetDEResults}
\title{
Obtain Differential Expression Analysis Results in a Two-condition Test 
}
\description{
Obtain DE analysis results in a two-condition test using the output of EBTest()
}
\usage{
GetDEResults(EBPrelim, FDR=0.05, Method="robust",
				     FDRMethod="hard", Threshold_FC=0.7,
					   Threshold_FCRatio=0.3, SmallNum=0.01)
}
\arguments{
	\item{EBPrelim}{Output from the function EBTest().}
	\item{FDR}{Target FDR, defaut is 0.05.}
	\item{FDRMethod}{"hard" or "soft".
	  Giving a target FDR alpha, either hard threshold and soft
	  threshold may be used.  If the hard threshold is preferred, DE transcripts are
	  defined as the the transcripts with PP(DE) greater than
	  (1-alpha). Using the hard threshold, any DE transcript in the list
	  has FDR <= alpha.
	
	  If the soft threshold is preferred, the DE transcripts are defined as the
	  transcripts with PP(DE) greater than crit_fun(PPEE, alpha). Using
	  the soft threshold, the list of DE transcripts has average FDR
	  alpha.
	
	  Based on results from our simulation studies, hard thresholds provide a better-controlled
	  empirical FDR when sample size is relatively small(Less than 10 samples in each condition).
	  User may consider the soft threshold when sample size is large to improve power.}
	\item{Method}{"robust" or "classic".
	  Using the "robust" option, EBSeq is more robust to genes with outliers and
	  genes with extremely small variances.
	  Using the "classic" option, the results will be more comparable to those obtained
	  by using the GetPPMat() function from earlier version (<= 1.7.0) of EBSeq.
	  Default is "robust".}
	\item{Threshold_FC}{Threshold for the fold change (FC) statistics.
	  The default is 0.7. The FC statistics are calculated as follows. 
		First the posterior FC estimates are calculated using PostFC() function.
	  The FC statistics is defined as exp(log(-|posterior FC|)) and therefore is always less than
		or equal to 1.
	  The default threshold was selected as the optimal threshold learned from our simulation studies. By setting the
	  threshold as 0.7, the expected FC for a DE transcript is less than 0.7
	  (or greater than 1/0.7=1.4).
	  User may specify their own threshold here. A higher (less conservative) threshold
	  may be used here when sample size is large. Our simulation results
	 indicated that when there are more than or equal to 5 samples in each condition,
	 a less conservative threshold will improve the power when the FDR is still well-controlled.
	        The parameter will be ignored if Method is set as "classic".}
	\item{Threshold_FCRatio}{Threshold for the fold change ratio (FCRatio) statistics.
	  The default is 0.3. The FCRatio statistics are calculated as follows. 
		First we get another revised fold change
	  statistic called Median-FC statistic for each transcript. 
		For each transcript, we calculate the median of
	  normalized expression values within each condition.
		The MedianFC is defined as exp(log(-|(C1Median+SmallNum)/(C2Median+SmallNum)|)).
		Note a small number is added to avoid Inf and NA. See SmallNum for more details.
		The FCRatio is calculated as exp(log(-|FCstatistics/MedianFC|)).
		Therefore it is always less than or equal to 1.
		The default threshold was selected as the optimal threshold learned from our simulation studies.	  
		By setting the threshold as 0.3, the FCRatio for a DE transcript is 
		expected to be larger than 0.3.
	}
	\item{SmallNum}{When calculating the FCRatio (or Median-FC), a small number is added for each transcript in each
	          condition to avoid Inf and NA. Default is 0.01.}
}
\details{
GetDEResults() function takes output from EBTest() function and output a list of 
DE transcripts under a target FDR. It also provides posterior probability estimates for each 
transcript.
}
\value{
  	\item{DEfound}{A list of DE transcripts.}
	  \item{PPMat}{Posterior probability matrix. Transcripts are following the same order as
		in the input matrix.
		Transcripts that were filtered by magnitude (in EBTest function), FC, or FCR
		are assigned with NA for both PPDE and PPEE.}
	  \item{Status}{Each transcript will be assigned with one of the following
		values: "DE", "EE", "Filtered: Low Expression", 
		"Filtered: Fold Change" and "Filtered: Fold Change Ratio".
		Transcripts are following the same order as in the input matrix.}
}
\references{
Ning Leng, John A. Dawson, James A. Thomson, Victor Ruotti, Anna I. Rissman, Bart M.G. Smits, Jill D. Haag, Michael N. Gould, Ron M. Stewart, and Christina Kendziorski. EBSeq: An empirical Bayes hierarchical model for inference in RNA-seq experiments. Bioinformatics (2013)
}
\author{
Ning Leng, Yuan Li
}
\seealso{
EBTest
}
\examples{
data(GeneMat)
str(GeneMat)
GeneMat.small = GeneMat[c(1:10,511:550),]
Sizes = MedianNorm(GeneMat.small)
EBOut = EBTest(Data = GeneMat.small,
	Conditions = as.factor(rep(c("C1","C2"), each = 5)),
	sizeFactors = Sizes, maxround = 5)
Out = GetDEResults(EBOut)
}
\keyword{ DE }
\keyword{ Two condition }
