\name{breastdata}
\alias{breastdata}
\docType{data}
\title{Breast cancer gene expression + DNA copy number data set from Chin et al
  (2006), Cancer Cell}
\description{
This data set consists of gene expression and DNA copy number
measurements on a set of 89 samples. This example is used in Witten,
Tibshirani and Hastie (2008).
}
\usage{data(breastdata)}
\format{
  The format is a list containing the following elements:
  - dna: a
  2149x89 matrix of CGH spots x Samples
  - rna: a 19672x89 matrix of
  Genes x Samples
  - chrom: a 2149-vector of chromosomal location of
  each CGH spot
  - nuc: a 2149-vector of nucleotide position for each
  CGH spot
  - gene: a 19672-vector wiith an accession number for each
  gene
  - genenames: a 19672-vector with a name for each gene
  - genechr: a 19672-vector with a chromosomal location for each
  gene
  - genedesc: a 19672-vector with a description for each gene
  - genepos: a 19672-vector with a nucleotide position for each gene
}
\details{
This data set can be used to perform integrative analysis of gene
expression and DNA copy number data, as in e.g. Witten and Tibshirani
and Hastie (2008). That is, we can look for sets of genes that are
associated with regions of chromosomal gain/loss.

Missing values were imputed using 5-nearest neighbors (see library pamr).
}
\source{
  This data set  was published in the following paper:

 Chin, K., DeVries, S., Fridlyand, J., Spellman, P., Roydasgupta, R.,
 Kuo,W.-L., Lapuk, A., Neve, R., Qian, Z., Ryder,
 T., Chen, F., Feiler, H., Tokuyasu, T., Kingsley, C., Dairkee, S.,
 Meng, Z., Chew, K., Pinkel, D., Jain, A., Ljung,
 B., Esserman, L., Albertson, D.,Waldman, F. & Gray, J. (2006), 'Genomic
 and transcriptional aberrations linked to breast cancer
 pathophysiologies', Cancer Cell 10, 529-541.

 It is publicly available at http://icbp.lbl.gov/breastcancer .
}
\references{
  Chin, DeVries, Fridlyand, et al. (2006) Cancer Cell 10, 529-541.

  Used as an example in Witten, Tibshirani and Hastie (2008) 'A
  penalized matrix decomposition, with applications to sparse principal
  components and canonical correlation analysis.' 
}
\examples{
data(breastdata)
attach(breastdata)
PlotCGH(dna[,1], chrom=chrom, main="Sample 1", nuc=nuc)
detach(breastdata)
}
\keyword{datasets}
