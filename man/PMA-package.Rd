\name{PMA-package}
\alias{PMA-package}
\alias{PMA}
\docType{package}
\title{
Penalized Multivariate Analysis
}
\description{
  This package is called PMA, for "Penalized Multivariate Analysis".
 It implements three methods:  A penalized matrix decomposition, sparse
 principal components analysis, and
sparse canonical correlations analysis. All are described in the paper "A penalized matrix
decomposition, with applications to sparse principal components and
canonical correlation analysis", by D Witten, R Tibshirani, and T
Hastie, and published in \emph{Biostatistics} (2009). 

The main functions are as follows:
(1) PMD
(2) CCA
(3) SPC

The first function, PMD, performs a penalized matrix decomposition.  CCA performs sparse canonical
correlation analysis. SPC
performs sparse principal components analysis.

There also are cross-validation functions for tuning parameter selection
for each of the above methods:
SPC.cv, PMD.cv,  CCA.permute. And PlotCGH results in nice plots for DNA copy number
data.

}
\details{
\tabular{ll}{
Package: \tab PMA\cr
Type: \tab Package\cr
Version: \tab 1.0.8\cr
Date: \tab 2011-08-11\cr
License: \tab GPL >= 2 \cr
}
}
\author{ Daniela M. Witten and Robert Tibshirani and Sam Gross and
  Balasubramanian Narasimhan
}
\references{ Witten, Tibshirani and Hastie (2009) A penalized matrix
  decomposition, with applications to sparse principal components and
  canonical correlation analysis. \emph{Biostatistics} 10(3): 515-534.
}

\keyword{ package }
