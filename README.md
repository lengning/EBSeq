# EBSeq
R/EBSeq is an R package for identifying genes and isoforms differentially expressed (DE) across 
two or more biological conditions in an RNA-seq experiment. Details can be found in Leng et al., 2013. 
The R/EBSeq package may be downloaded below. A vignette is also available there. 
It provides the syntax required for identifying DE genes and isoforms in a two-group RNA-seq experiment 
as well for identifying DE genes across more than two conditions (as noted in the vignette, the commands 
for identifying DE isoforms across more than two conditions are the same as those required for gene-level analysis).


Installation

EBSeq is currently available at Bioconductor (EBSeq Bioconductor page). To install, type the following commands in R:

source("http://bioconductor.org/biocLite.R")

devel = "http://bioconductor.org/packages/3.0/bioc"

biocLite("EBSeq", siteRepos = devel, type="source")

library(EBSeq)

packageVersion("EBSeq")
