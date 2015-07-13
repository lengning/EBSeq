# EBSeq Q & A


## ReadIn data

csv file:

```
In=read.csv("FileName", stringsAsFactors=F, row.names=1, header=T)

Data=data.matrix(In)
```

txt file:
```
In=read.table("FileName", stringsAsFactors=F, row.names=1, header=T)

Data=data.matrix(In)
```
check str(Data) and make sure it is a matrix instead of data frame. You may need to play around with the row.names and header option depends on how the input file was generated.



## GetDEResults() function not found

You may on an earlier version of EBSeq. The GetDEResults function
was introduced since version 1.7.1.
The latest release version could be found at:
http://www.bioconductor.org/packages/devel/bioc/html/EBSeq.html
And you may check your package version by typing packageVersion("EBSeq")


## Visualizing DE genes/isoforms

To generate a heatmap, you may consider the heatmap.2 function in gplots package.
For example, you may run
```
heatmap.2(NormalizedMatrix[GenesOfInterest,], scale="row", trace="none", Colv=F)
```
The normalized matrix may be obtained from GetNormalizedMat() function.


## My favorite gene/isoform has NA in PP (status "NoTest")

The NoTest status comes from two sources

1) Using the default parameter settings of EBMultiTest(), the function
will not test on genes with more than 75% values < 10 to ensure better
model fitting. To disable this filter, you may set Qtrm=1 and
QtrmCut=0.

2) numerical over/underflow in R. That happens when the within
condition variance is extremely large or small. I did implemented a numerical
approximation step to calculate the approximated PP for these genes
with over/underflow. Here I use 10^-10 to approximate the parameter p
in the NB distribution for these genes (I set it to a small value
since I want to cover more over/underflow genes with low
within-condition variation). You may try to tune this value (to a larger value) in the
approximation by setting ApproxVal in EBTest() or EBMultiTest() function. 

## Can I run more than 5 iterations when running EBSeq via RSEM wrapper?

Yes you may modify the script rsem-for-ebseq-find-DE under RSEM/EBSeq
change line 36
```
EBOut <- EBTest(Data = DataMat, NgVector = ngvector, Conditions =
conditions, sizeFactors = Sizes, maxround = 5)
```
to
```
EBOut <- EBTest(Data = DataMat, NgVector = ngvector, Conditions =
conditions, sizeFactors = Sizes, maxround = 10)
```
If you are running multiple condition analysis, you will need to change line 53:
```
MultiOut <- EBMultiTest(Data = DataMat, NgVector = ngvector,
Conditions = conditions, AllParti = patterns, sizeFactors = Sizes,
maxround = 5)
```
```
MultiOut <- EBMultiTest(Data = DataMat, NgVector = ngvector,
Conditions = conditions, AllParti = patterns, sizeFactors = Sizes,
maxround = 10)
```
You will need to redo make after you make the changes.

## I saw a gene has significant FC but is not called as DE by EBSeq, why does that happen?

EBSeq calls a gene as DE (assign high PPDE) if the across-condition variability is significantly larger than the within-condition
variability. In the cases that a gene has large within-condition variation, although the FC across two conditions is large (small), 
the across-condition difference could still be explained by biological variation within condition. In these cases the gene/isoform
will have a moderate PPDE.
