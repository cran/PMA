\name{CCA}
\alias{CCA}
\title{Perform sparse canonical correlation analysis using the penalized
matrix decomposition.}
\description{
Given matrices X and Z, which represent two sets of features on the same
set of samples, find sparse u and v such that u'X'Zv is large. If the
columns of Z are ordered (and type="ordered") then v will also be
smooth. For X and Z, the
samples are on the rows and the features are on the columns. X and Z
must have same number of rows, but may (and usually will) have different
numbers of columns. 
} 
\usage{
CCA(x, z, type=c("standard", "ordered"), sumabs=.4, sumabsu=4,
sumabsv=NULL,
lambda=NULL, K=1, niter=25,v=NULL, trace=TRUE, standardize=TRUE, xnames=NULL,
znames=NULL, chrom=NULL, upos=FALSE, uneg=FALSE, vpos=FALSE, vneg=FALSE,
outcome=NULL,y=NULL,cens=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Data matrix; samples are rows and columns are
    features. Cannot contain missing values.}
  \item{z}{Data matrix; samples are rows and columns are
      features. Note that x and z must have the same number of rows, but
      may (and generally will) have different numbers of columns. Cannot
      contain missing values.}
    \item{type}{Are the columns of z unordered (type="standard") or
      ordered (type="ordered")? If "standard", then a lasso penalty is
      applied to v, to enforce sparsity. If "ordered" (generally used
      for CGH data), then a fused
      lasso penalty is applied, to enforce both sparsity and smoothness.}
   \item{sumabs}{Used only for type "standard". Controls both sumabsu and sumabsv at once. Must be
     between 0 and 1. If either sumabsu
     or sumabsv is NULL, then the value of sumabs will be used to
     determine values for sumabsu and sumabsv, as follows: sumabsu will
     be set to sqrt(ncol(x))*sumabs and sumabsv will be set to
     sqrt(ncol(z))*sumabs. Default value is 1.5. If sumabsu and sumabsv
     both are non-null, then this argument is ignored.}
  \item{sumabsu}{How sparse do you want u to be? This is the sum of
    absolute values of elements of u. It must be between 1 and the
    square root of the number of columns of x. The smaller it is,
  the sparser u will be. Default is NULL (sumabs is used instead).}
  \item{sumabsv}{Used only for type "standard". How sparse do you want v to be? This is the sum of
    absolute values of elements of v. It must be between 1 and square
    root of number of columns of z. The smaller it is, the sparser v
    will be. Default is NULL (sumabs is used instead).}
  \item{lambda}{Used only for type "ordered", controls fused lasso
    penalty on v, which takes the form $lambda( sum_j |v_j| + sum_j
    |v_j - v_(j-1)|)$.}
    \item{K}{How many factors u and v do you want? Default is 1.}
  \item{niter}{How many iterations should be performed? Default is 25.}
  \item{v}{The first K columns of the v matrix of the SVD of X'Z. If
    NULL, then the SVD of X'Z will be computed inside the CCA function. However, if
    you plan to run this function multiple times, then save a copy of
    this argument so that it does not need to be re-computed (since that
    process can be time-consuming if X and Z both have high dimension).}
  \item{trace}{Print out progress?}
  \item{standardize}{Should the columns of x and z be centered (to have mean zero)
    and scaled (to have standard deviation 1)? Default is TRUE.}
  \item{xnames}{An optional vector of column names for x.}
  \item{znames}{An optional vector of column names for z.}
  \item{chrom}{Used only if type is "ordered"; allows user to specify a
    vector of length ncol(z) giving the chromosomal location of each CGH
    spot. This is so that smoothness will be enforced within each
    chromosome, but not between chromosomes.}
  \item{upos}{If TRUE, then require elements of u to be positive. FALSE
    by default.}
  \item{uneg}{If TRUE, then require elements of u to be negative. FALSE
    by default.}
  \item{vpos}{Only applicable if type="standard". If TRUE, require
    elements of v to be positive. FALSE by default.}
  \item{vneg}{Only applicable if type="standard". If TRUE, require
    elements of v to be negative. FALSE by default.}
  \item{outcome}{If you would like to incorporate a phenotype into CCA
    analysis - that is, you wish to find features that are correlated
    across the two data sets and also correlated
    with a phenotype - then use one of "survival", "multiclass", or
    "quantitative" to indicate outcome type. Default is NULL.}
  \item{y}{If outcome is not NULL, then this is a vector of phenotypes -
    one for each row of x and z. If outcome is "survival" then these are
    survival times; must be non-negative. If outcome is "multiclass"
    then these are class labels (1,2,3,...). Default NULL.}
  \item{cens}{If outcome is "survival" then these are censoring statuses
    for each observation. 1 is complete, 0 is censored. Default NULL.}
}
\details{
This function is useful for performing an integrative analysis of two
data sets taken on the same set of samples: for instance, gene
expression and CGH measurements on the same set of patients. It takes in
two data sets, called x and z, each of which have (the same set of)
samples on the rows. If z is a matrix of CGH data with *ordered* CGH
spots on the columns, then use type="ordered". If z consists   of
unordered columns, then use type="standard".
  
  This function performs the penalized matrix decomposition on the data
  matrix $X'Z$. Therefore, the results should be the same  as running
  the PMD function on t(x)\%*\%z. However, when ncol(x)>>nrow(x) and
  ncol(z)>>nrow(z) then using the CCA function is much faster because it
  avoids computation of $X'Z$. (It is easy to verify that
  PMD(t(x)\%*\%z, sumabs=.5, type="standard", center=FALSE) gives same results as
  CCA(x, z, sumabs=.5, type="standard", standardize=FALSE) for some
  matrices x and z).

  The CCA criterion is as follows: find unit vectors $u$ and $v$ such
  that $u'X'Zv$ is maximized subject to constraints on $u$ and $v$.   If
  type="standard" then the constraints on $u$ and $v$ are lasso
  ($L_1$). If type="ordered" then the constraint on $u$ is a lasso
  constraint, and there is a fused lasso constraint on $v$ (promoting
  sparsity and smoothness).
  

  When type is "standard": If either sumabsu or sumabsv is NULL, then sumabs must be non-NULL. In
  this case, sumabs will be used to set values for both sumabsu and
  sumabsv, as follows: sumabsu will be set to sumabs*sqrt(ncol(x)) and
  sumabsv will be set to sumabs*sqrt(ncol(z)). The sumabs argument is
  ignored if sumabsu and sumabsv both are non-NULL.

  When type is "ordered": lambda controls the amount of sparsity and
  smoothness in v, via the fused lasso penalty: $lambda sum_j |v_j| +
  lambda sum_j |v_j - v_(j-1)|$. If NULL, then it will be chosen
  adaptively from the data. Sumabsu is the bound on the sum of absolute values
  of elements of u - the greater it is, the sparser u will be.

  When running CCA on gene expression and CGH data (type="ordered"),
  there are multiple ways to do it:

  (1) using all gene expression data and all CGH data, regardless of
  chromosome.
  
  (2) using CGH data on chromosome j and gene expression data on all
  chromosomes. This will allow for discovery of cis/trans interactions
  relating to copy number change on chrom j. Then, repeat for other
  chromosomes.
  
  (3) using CGH data on chromosome j and gene expression data on
  chromosome j. Cis interactions can be found this way.
  
  (4) using CGH data on chromosome j and gene expression data NOT on
  chromosome j; one can find trans interactions this way...
}
\value{
  \item{u}{u is output. If you asked for multiple factors then each
    column of u is a factor. u has dimension nxK if you asked for K factors.}
  \item{v}{v is output. If you asked for multiple factors then each
    column of v is a factor. v has dimension pxK if you asked for K
    factors.}
  \item{d}{A vector of length K, which can alternatively be computed as
    the diagonal of the matrix $u'X'Zv$.}
  \item{v.init}{The first K factors of the v matrix of the SVD of
    x'z. This is saved in case this function will be re-run later.}
}
\references{Witten, DM and Tibshirani, R and T Hastie (2008) A penalized
  matrix decomposition, with applications to
  sparse principal components and canonical correlation
  analysis. Submitted. <http://www-stat.stanford.edu/~dwitten>}
\author{Daniela M. Witten and Robert Tibshirani}
\seealso{\link{PMD},\link{CCA.permute}}
\examples{
# first, do CCA with type="standard"
# A simple simulated example
u <- matrix(c(rep(1,25),rep(0,75)),ncol=1)
v1 <- matrix(c(rep(1,50),rep(0,450)),ncol=1)
v2 <- matrix(c(rep(0,50),rep(1,50),rep(0,900)),ncol=1)
x <- u\%*\%t(v1) + matrix(rnorm(100*500),ncol=500)
z <- u\%*\%t(v2) + matrix(rnorm(100*1000),ncol=1000)
# Can run CCA with default settings, and can get e.g. 3 components
out <- CCA(x,z,type="standard",K=3)
print(out,verbose=TRUE) # To get less output, just print(out)
# Or can use CCA.permute to choose optimal parameter values
perm.out <- CCA.permute(x,z,type="standard",nperms=7,sumabss=seq(0.1,.75,len=12))
print(perm.out)
plot(perm.out)
out <- CCA(x,z,type="standard",K=1,sumabs=perm.out$bestsumabs, v=perm.out$v.init)
print(out)
# Now try CCA with a constraint that elements of u must be negative and
# elements of v must be positive:
perm.out <- CCA.permute(x,z,type="standard",nperms=7,
sumabss=seq(0.1,.75,len=12), uneg=TRUE, vpos=TRUE)
print(perm.out)
plot(perm.out)
out <- CCA(x,z,type="standard",K=1,sumabs=perm.out$bestsumabs,
v=perm.out$v.init, uneg=TRUE, vpos=TRUE)
print(out)

# Suppose we also have a quantitative outcome, y, and we want to find
# features in x and z that are correlated with each other and with the
# outcome:
y <- rnorm(nrow(x))
perm.out <- CCA.permute(x,z,type="standard",outcome="quantitative",y=y)
print(perm.out)
out<-CCA(x,z,type="standard",outcome="quantitative",y=y,sumabs=perm.out$bestsumabs)
print(out)


# now, do CCA with type="ordered"
# Example involving the breast cancer data: gene expression + CGH
set.seed(22)
data(breastdata)
attach(breastdata)
dna <- t(dna)
rna <- t(rna)
perm.out <- CCA.permute(rna,dna[,chrom==1],type="ordered",nperms=5)
# We run CCA using all gene exp. data, but CGH data on chrom 1 only.
print(perm.out)
plot(perm.out)
out <- CCA(rna,dna[,chrom==1], type="ordered",sumabsu=perm.out$bestsumabsu,
v=perm.out$v.init, lambda=perm.out$lambda, xnames=substr(genedesc,1,20),
znames=paste("Pos", sep="", nuc[chrom==1])) # Save time by inputting  lambda and v
print(out) # could do print(out,verbose=TRUE)
print(genechr[out$u!=0]) # Cool! The genes associated w/ gain or loss
# on chrom 1 are located on chrom 1!!
par(mfrow=c(1,1))
PlotCGH(out$v, nuc=nuc[chrom==1], chrom=chrom[chrom==1],
main="Regions of gain/loss on Chrom 1 assoc'd with gene expression")
detach(breastdata)

}

