\name{PlotCGH}
\alias{PlotCGH}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot CGH data}
\description{
Given a vector of gains/losses at CGH spots, this makes a plot of
gain/loss on each chromosome.
}
\usage{
PlotCGH(array,chrom=NULL,nuc=NULL,main="",scaleEachChrom=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{array}{A vector containing the chromosomal location of each CGH spot.}
 \item{chrom}{A numeric vector of the same length as "array"; its values should
   indicate the
   chromosome that each CGH spot is on (for instance, for human genomic
   data, values of chrom should range
   from 1 to 24). If NULL, then it is assumed that all elements of
   'array' are on the same chromosome.}
 \item{nuc}{A numeric vector of same length as "array", indicating the
   nucleotide position of each CGH spot. If NULL, then the function
   assumes that each CGH spot corresponds to a consecutive position. E.g. if
   there are 200 CGH spots on chromosome 1, then they are located at
   positions 1,2,...,199,200.}
 \item{main}{Give your plot a title.}
 \item{scaleEachChrom}{Default is TRUE. This means that each chromosomes
 CGH spots are divided by 1.1 times the max of the CGH spots on that
 chromosome. This way, the CGH spots on each chromosome of the plot are
 as big as possible (i.e. easy to see). If FALSE, then all of the CGH
 spots are divided by 1.1 times the max of ALL the CGH spots. This means
that on some chromosomes CGH spots might be hard to see, but has the
advantage that now relative magnitudes of CGH spots on different chromosomes can
be seen from figure.}
}
\details{
This function makes a plot of regions of genomic gain/loss. 
}
\references{Witten DM, Tibshirani R and T Hastie (2008) A penalized matrix
  decomposition with applications to sparse principal components and
  canonical correlation analysis. Submitted. <http://www-stat.stanford.edu/~dwitten>}
\author{Daniela M. Witten (adapted from Pei Wang and Rob Tibshirani's
  cghFLasso package)}
\seealso{\link{PMD}, \link{PMD.cv}, \link{CCA}, \link{CCA.permute}}
\examples{
# Use breast data
data(breastdata)
attach(breastdata)

# dna contains CGH data and chrom contains chromosome of each CGH spot;
# nuc contains position of each CGH spot.
dna <- t(dna)
PlotCGH(dna[1,],chrom=chrom,nuc=nuc,main="Sample 1: All Chromosomes")
PlotCGH(dna[1,chrom==1], chrom=chrom[chrom==1], nuc=nuc[chrom==1],
main= "Sample 1: Chrom 1")
PlotCGH(dna[1,chrom<=3], chrom=chrom[chrom<=3], nuc=nuc[chrom<=3],
 main="Sample 1: Chroms 1, 2, and 3")
detach(breastdata)
}

