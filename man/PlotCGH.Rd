% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PlotFusedLasso.R
\name{PlotCGH}
\alias{PlotCGH}
\title{Plot CGH data}
\usage{
PlotCGH(array, chrom = NULL, nuc = NULL, main = "", scaleEachChrom = TRUE)
}
\arguments{
\item{array}{A vector containing the chromosomal location of each CGH spot.}

\item{chrom}{A numeric vector of the same length as "array"; its values
should indicate the chromosome that each CGH spot is on (for instance, for
human genomic data, values of chrom should range from 1 to 24). If NULL,
then it is assumed that all elements of 'array' are on the same chromosome.}

\item{nuc}{A numeric vector of same length as "array", indicating the
nucleotide position of each CGH spot. If NULL, then the function assumes
that each CGH spot corresponds to a consecutive position. E.g. if there are
200 CGH spots on chromosome 1, then they are located at positions
1,2,...,199,200.}

\item{main}{Give your plot a title.}

\item{scaleEachChrom}{Default is TRUE. This means that each chromosomes CGH
spots are divided by 1.1 times the max of the CGH spots on that chromosome.
This way, the CGH spots on each chromosome of the plot are as big as
possible (i.e. easy to see). If FALSE, then all of the CGH spots are divided
by 1.1 times the max of ALL the CGH spots. This means that on some
chromosomes CGH spots might be hard to see, but has the advantage that now
relative magnitudes of CGH spots on different chromosomes can be seen from
figure.}
}
\description{
Given a vector of gains/losses at CGH spots, this makes a plot of gain/loss
on each chromosome.
}
\details{
This function makes a plot of regions of genomic gain/loss.
}
\examples{
\dontrun{
# Use breast data
breastdata <- download_breast_data()
with(breastdata, {
# dna contains CGH data and chrom contains chromosome of each CGH spot;
# nuc contains position of each CGH spot.
dna <- t(dna)
ch1 <- which(chrom == 1)
PlotCGH(dna[1,],chrom=chrom,nuc=nuc,main="Sample 1: All Chromosomes")
PlotCGH(dna[1,ch1], chrom=chrom[ch1], nuc=nuc[ch1],
main= "Sample 1: Chrom 1")
chlt3 = which(chrom <= 3)
PlotCGH(dna[1,chlt3], chrom=chrom[chlt3], nuc=nuc[chlt3],
 main="Sample 1: Chroms 1, 2, and 3")
} )
}
}
\references{
Witten D. M., Tibshirani R.,  and Hastie, T. (2009)
\emph{A penalized matrix decomposition, with applications to sparse principal components and canonical correlation analysis}, \emph{Biostatistics, Gol 10 (3), 515-534, Jul 2009}\cr
}
\seealso{
\link{PMD}, \link{PMD.cv}, \link{CCA}, \link{CCA.permute}
}
