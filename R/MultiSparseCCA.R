fastsvd <- function(x,z){
  # fast svd of t(x)%*%z, where ncol(x)>>nrow(x) and same for z
  xx=x%*%t(x)
  xx2=msqrt(xx)
  y=t(z)%*%xx2
  a=svd(y)
  v=a$u
  d=a$d
  zz=z%*%t(z)
  zz2=msqrt(zz)
  y=t(x)%*%zz2
  a=svd(y)
  u=a$u
  return(list(u=u,v=v,d=d))
}

msqrt <- function(x){
  eigenx <- eigen(x)
  return(eigenx$vectors%*%diag(sqrt(pmax(0,eigenx$values)))%*%t(eigenx$vectors))
}




SparseCCA <- function(x,y,v,sumabsu,sumabsv,niter,trace=TRUE, upos, uneg, vpos){
  if(sumabsu<1 || sumabsv<1 || sumabsu>sqrt(ncol(x)) || sumabsv>sqrt(ncol(y))) stop("Sumabsu and Sumabsv must be between 1 and sqrt(ncol(x)), sqrt(ncol(z))")
  vold <- rnorm(length(v))
  for(i in 1:niter){
    if(sum(abs(vold-v))>1e-6){
      if(trace) cat(i,fill=F)
      # Update u #
      argu <- t(x)%*%(y%*%v)
      if(upos) argu <- pmax(argu,0)
      if(uneg) argu <- pmin(argu,0)
      lamu <- BinarySearch(argu,sumabsu)
      su <- soft(argu,lamu)
      u <-  matrix(su/l2n(su), ncol=1)
      # Done updating u #
      # Update v #
      vold <- v
      argv <- (t(u)%*%t(x))%*%y
      if(vpos) argv <- pmax(argv,0)
      lamv <- BinarySearch(argv,sumabsv)
      sv <- soft(argv, lamv)
      v <-  matrix(sv/l2n(sv),ncol=1)
      # Done updating v #
    }
  }
  if(trace) cat(fill=T)
  # Update d #
  d <-  sum((x%*%u)*(y%*%v))
  # Done updating d #
  if(sum(is.na(u))>0 || sum(is.na(v))>0){
    u <- matrix(rep(0,ncol(x)),ncol=1)
    v <- matrix(rep(0,ncol(y)),ncol=1)
    d <- 0
  }
  return(list(u=u,v=v,d=d))
}

CheckVs <- function(v,x,z,K){ # If v is NULL, then get v as apropriate.
  if(!is.null(v) && !is.matrix(v)) v <- matrix(v,nrow=ncol(z))
  if(!is.null(v) && ncol(v)<K) v <- NULL
  if(!is.null(v) && ncol(v)>K) v <- matrix(v[,1:K],ncol=K)
  if(is.null(v) && ncol(z)>nrow(z) && ncol(x)>nrow(x)){
    v <- matrix(fastsvd(x,z)$v[,1:K],ncol=K)
  } else if (is.null(v) && (ncol(z)<=nrow(z) || ncol(x)<=nrow(x))){
    v <- matrix(svd(t(x)%*%z)$v[,1:K],ncol=K)
  }
  return(v)
}


ftrans <- function(a){
  return(log((1+a)/(1-a)))
}

CCAL1L1.permute <- function(x,z,sumabss=seq(.2,.6,len=10),niter=3,v=NULL,trace=TRUE,nperms=25, standardize=TRUE, upos, uneg, vpos){
  call <- match.call()
  if(standardize){
    x <- scale(x,TRUE,TRUE)
    z <- scale(z,TRUE,TRUE)
  } 
  v <- CheckVs(v,x,z,1)
  ccperms=nnonzerous.perms=nnonzerovs.perms=matrix(NA, length(sumabss), nperms)
  ccs=nnonzerous=nnonzerovs=numeric(length(sumabss))
  for(i in 1:nperms){
    if(trace) cat("\n Permutation ",i," out of ", nperms, " ")
    samp <- sample(1:nrow(z))
    for(sabs in 1:length(sumabss)){
      if(trace) cat(sabs,fill=FALSE)
      if(i==1){
        out <- CCAL1L1(x,z,sumabs=sumabss[sabs],niter=niter,v=v,trace=FALSE, upos=upos, uneg=uneg, vpos=vpos)
        nnonzerous[sabs] <- sum(out$u!=0)
        nnonzerovs[sabs] <- sum(out$v!=0)
        if(mean(out$u==0)!=1 && mean(out$v==0)!=1) ccs[sabs] <- cor(x%*%out$u,z%*%out$v)
      }
      out <- CCAL1L1(x,z[samp,],sumabs=sumabss[sabs],niter=niter,v=v,trace=FALSE, upos=upos, uneg=uneg, vpos=vpos)
      nnonzerous.perms[sabs,i] <- sum(out$u!=0)
      nnonzerovs.perms[sabs,i] <- sum(out$v!=0)
      if(mean(out$u==0)!=1 && mean(out$v==0)!=1){
        ccperms[sabs,i] <- cor(x%*%out$u,z[samp,]%*%out$v)
      } else {
        ccperms[sabs,i] <- 0
      }
    }
  }
  cc.norm <- ftrans(ccs)
  ccperm.norm <- ftrans(ccperms)
  zstats <- (cc.norm - rowMeans(ccperm.norm))/(apply(ccperm.norm,1,sd) + .05)
    # 0.05 added to the denominator to avoid getting zstat of INFINITY
  if(trace) cat(fill=T)
  pvals <- apply(sweep(ccperms,1,ccs,"-")>=0,1,mean)
  results <- list(zstats=zstats,sumabss=sumabss,bestsumabss=sumabss[which.max(zstats)], cors=ccs, corperms=ccperms, ft.cors=cc.norm,ft.corperms=rowMeans(ccperm.norm),nnonzerous=nnonzerous,nnonzerovs=nnonzerovs, nnonzerous.perm=rowMeans(nnonzerous.perms),nnonzerovs.perm=rowMeans(nnonzerovs.perms),call=call,v.init=v,pvals=pvals,nperms=nperms)
  class(results) <- "ccaperml1l1"
  return(results)
}

plot.ccaperml1l1 <- function(x,...){
  sumabss <- x$sumabss
  ccs <- x$cors
  nperms <- x$nperms
  zstats <- x$zstats
  ccperms <- x$corperms
  par(mfrow=c(2,1))
  plot(sumabss, ccs, main="Correlations For Real/Permuted Data", xlab="Value of sumabs", ylab="Correlations", ylim=range(ccperms,ccs))
  points(sumabss,ccs,type="l")
  for(i in 1:nperms){
    points(sumabss,ccperms[,i],col="green")
  }
  plot(sumabss,zstats,main="Z", xlab="Value of sumabs", ylab="Z score")
  lines(sumabss,zstats)
}

CCAL1FL.permute <- function(x,z,sumabsus=seq(2,8,len=10),lambda=NULL,niter=3,v=NULL,trace=TRUE,nperms=25, standardize=TRUE,chrom=NULL,nuc=NULL, upos=FALSE, uneg=FALSE){
  call <- match.call()
  if(standardize){
    x <- scale(x,TRUE,TRUE)
    z <- scale(z, TRUE, TRUE)
  }
  v <- CheckVs(v,x,z,1)
  if(is.null(lambda)) lambda <- ChooseLambda1Lambda2(as.numeric(v))
  ccperms=nnonzerous.perms=nnonzerovs.perms=matrix(NA, length(sumabsus), nperms)
  ccs=nnonzerous=nnonzerovs=numeric(length(sumabsus))
  storevs <- NULL
  for(i in 1:nperms){
    if(trace) cat("\n Permutation ",i," out of ", nperms, " ")
    samp <- sample(1:nrow(z))
    for(sabs in 1:length(sumabsus)){
      if(trace) cat(sabs,fill=FALSE)
      if(i==1){
        out <- CCAL1FL(x,z,sumabsu=sumabsus[sabs],lambda=lambda,niter=niter,v=v,trace=FALSE,K=1,chrom=chrom, upos=upos, uneg=uneg)
        nnonzerous[sabs] <- sum(out$u!=0)
        nnonzerovs[sabs] <- sum(out$v!=0)
        if(mean(out$u==0)!=1 && mean(out$v==0)!=1) ccs[sabs] <- cor(x%*%out$u,z%*%out$v)
        storevs <- cbind(storevs, out$v)
      }
      out <- CCAL1FL(x,z[samp,],sumabsu=sumabsus[sabs],lambda=lambda,niter=niter,v=v,trace=FALSE,K=1,chrom=chrom, upos=upos, uneg=uneg)
      nnonzerous.perms[sabs,i] <- sum(out$u!=0)
      nnonzerovs.perms[sabs,i] <- sum(out$v!=0)
      if(mean(out$u==0)!=1 && mean(out$v==0)!=1){
        ccperms[sabs,i] <- cor(x%*%out$u,z[samp,]%*%out$v)
      } else {
        ccperms[sabs,i] <- 0
      }
    }
  }
  cc.norm <- ftrans(ccs)
  ccperm.norm <- ftrans(ccperms)
  zstats <- (cc.norm - rowMeans(ccperm.norm))/(apply(ccperm.norm,1,sd) + .05)
  if(trace) cat(fill=T)
  pvals <- apply(sweep(ccperms,1,ccs,"-")>=0,1,mean)
  results <- list(zstats=zstats,sumabsus=sumabsus,bestsumabsu=sumabsus[which.max(zstats)], lambda=lambda,cors=ccs, corperms=ccperms, ft.cors=cc.norm,ft.corperms=rowMeans(ccperm.norm),nnonzerous=nnonzerous,nnonzerovs=nnonzerovs, nnonzerous.perm=rowMeans(nnonzerous.perms),nnonzerovs.perm=rowMeans(nnonzerovs.perms),call=call,v.init=v, pvals=pvals,nperms=nperms,chrom=chrom,nuc=nuc,storevs=storevs)
  class(results) <- "ccaperml1fl"
  return(results)
}

plot.ccaperml1fl <- function(x,...){
  sumabsus <- x$sumabsus
  nperms <- x$nperms
  zstats <- x$zstats
  chrom <- x$chrom
  nuc <- x$nuc
  ccs <- x$cors
  ccperms <- x$corperms
  storevs <- x$storevs
  par(mfrow=c(3,1))
  plot(sumabsus, ccs, main="Correlations For Real/Permuted Data", xlab="Value of sumabsu", ylab="Correlations", ylim=range(ccperms,ccs))
  points(sumabsus,ccs,type="l")
  for(i in 1:nperms){
    points(sumabsus,ccperms[,i],col="green")
  }
  plot(sumabsus,zstats,main="Z", xlab="Value of sumabsu", ylab="Z score")
  lines(sumabsus,zstats)
  if(is.null(chrom)) chrom <- rep(1, nrow(storevs))
  PlotCGH(storevs[,which.max(zstats)], chrom=chrom, nuc=nuc, main="Best V obtained")
}

print.ccaperml1l1 <- function(x,...){
  cat("Call: ")
  dput(x$call)
  cat("\n\n")
  tab <- round(cbind(round(x$sumabss,3), x$zstats,x$pvals,x$cors,rowMeans(x$corperms),x$ft.cors,x$ft.corperms,x$nnonzerous,x$nnonzerovs),3)
  dimnames(tab) <- list(1:length(x$sumabss), c("Sumabss", "Z","P-Value","Cors","Cors Perm", "FT(Cors)", "FT(Cors Perm)", "# U's non-zero", "# V's non-zero"))
  print(tab,quote=F)
  cat("Highest z score: ", max(x$zstats), "\n")
  cat("Sumabs corresponding to highest z score: ", x$bestsumabs, "\n")
  if(x$upos) cat("U's required to be positive.", fill=TRUE)
  if(x$uneg) cat("U's required to be negative.", fill=TRUE)
  if(x$vpos) cat("V's required to be positive.", fill=TRUE)
}

print.ccaperml1fl <- function(x,...){
  cat("Call: ")
  dput(x$call)
  cat("\n\n")
  tab <- round(cbind(round(x$sumabsus,3), x$zstats,x$pvals,x$cors,rowMeans(x$corperms),x$ft.cors,x$ft.corperms,x$nnonzerous,x$nnonzerovs),3)
  dimnames(tab) <- list(1:length(x$sumabsus), c("Sumabsus", "Z","P-Value","Cors","Cors Perm", "FT(Cors)", "FT(Cors Perm)", "# U's non-zero", "# V's non-zero"))
  print(tab,quote=F)
  cat("Highest z score: ", max(x$zstats), "\n")
  cat("Sumabsu corresponding to highest z score: ", x$bestsumabsu, "\n")
  cat("Lambda used: ", x$lambda, "\n")
  if(x$upos) cat("U's required to be positive.", fill=TRUE)
  if(x$uneg) cat("U's required to be negative.", fill=TRUE)
}

CCA <- function(x, z, type=c("standard", "ordered"), sumabs=.5, sumabsu=4, sumabsv=NULL, lambda=NULL, K=1, niter=25, v=NULL, trace=TRUE, standardize=TRUE, xnames=NULL, znames=NULL, chrom=NULL, upos=FALSE, uneg=FALSE, vpos=FALSE){
  if(upos && uneg) stop("At most one of upos and uneg should be TRUE!")
  if(type=="ordered" && vpos) stop("Cannot require elements of v to be positive if type is ordered")
  call <- match.call()
  type <- match.arg(type)
  if(type=="standard"){
    out <- CCAL1L1(x,z,sumabs=sumabs,sumabsu=sumabsu, sumabsv=sumabsv,K=K,niter=niter,v=v,trace=trace,standardize=standardize,xnames=xnames,znames=znames,upos=upos,uneg=uneg,vpos=vpos)
    class(out) <- "ccal1l1"
  } else if (type=="ordered"){
    out <- CCAL1FL(x,z,K=K,sumabsu=sumabsu,lambda=lambda,chrom=chrom,niter=niter,v=v,trace=trace,standardize=standardize,xnames=xnames,znames=znames,upos=upos,uneg=uneg)
    class(out) <- "ccal1fl"
  }
  out$call <- call
  return(out)
}

CCA.permute <- function(x,z,type=c("standard", "ordered"), sumabss=seq(.1,.9,len=10),sumabsus=NULL, lambda=NULL,niter=3,v=NULL,trace=TRUE,nperms=25, standardize=TRUE, chrom=NULL, nuc=NULL, upos=FALSE, uneg=FALSE, vpos=FALSE){
  if(type=="ordered" && vpos) stop("If type=ordered then you cannot require elements of v to be positive!")
  if(is.null(sumabsus) && type=="ordered") sumabsus <- pmax(1.2, exp(seq(log(.01*sqrt(ncol(x))), log(.9*sqrt(ncol(x))), len=10)))
  call <- match.call()
  type <- match.arg(type)
  if(type=="standard"){
    out <- CCAL1L1.permute(x=x,z=z,sumabss=sumabss,niter=niter, v=v, trace=trace, nperms=nperms, standardize=standardize, upos=upos, uneg=uneg, vpos=vpos)
    class(out) <- "ccaperml1l1"
  } else if (type=="ordered"){
    out <- CCAL1FL.permute(x=x,z=z,sumabsus=sumabsus, lambda=lambda,niter=niter,v=v,trace=trace,nperms=nperms,standardize=standardize,chrom=chrom,nuc=nuc, upos=upos, uneg=uneg)
    class(out) <- "ccaperml1fl"
  }
  out$call <- call
  out$upos <- upos
  out$uneg <- uneg
  out$vpos <- vpos
  return(out)
}

CCAL1L1 <- function(x,z,sumabs=.5,sumabsu=NULL,sumabsv=NULL,K=1,niter=25,v=NULL, trace=TRUE, standardize=TRUE, xnames=NULL, znames=NULL,upos=FALSE,uneg=FALSE,vpos=FALSE){
  call <- match.call()
  if(sum(is.na(x))+sum(is.na(z)) > 0) stop("Cannot have NAs in x or z")
  if(nrow(x)!=nrow(z)) stop("x and z must have same number of rows")
  if(standardize){
    x <- scale(x,TRUE,TRUE)
    z <- scale(z,TRUE,TRUE)
  }
  v <- CheckVs(v,x,z,K)
  if(is.null(sumabs)&&(is.null(sumabsu)||is.null(sumabsv))) stop("must enter sumabs, or sumabsu&sumabsv.")
  if(is.null(sumabsu)||is.null(sumabsv)){
#    if(!is.null(sumabsu) || !is.null(sumabsv)) warning("Since sumabsu and sumabsv weren't both given, both are being ignored and sumabs is being used instead.")
    sumabsu <- sumabs*sqrt(ncol(x))
    sumabsv <- sumabs*sqrt(ncol(z))
  }
  out <- MultiSparseCCA(x,z,sumabsu=sumabsu,sumabsv=sumabsv,niter=niter,K=K,v=v,trace=trace,upos=upos,uneg=uneg,vpos=vpos)
  out$call <- call
  out$xnames <- xnames
  out$znames <- znames
  class(out) <- "ccal1l1"
  return(out)
}

print.ccal1l1 <- function(x,verbose=FALSE,...){
  cat("Call: ")
  dput(x$call)
  cat("\n\n")
  cat("Num non-zeros u's: ", apply(x$u!=0,2,sum), "\n")
  cat("Num non-zeros v's: ", apply(x$v!=0,2,sum), "\n")
  cat("Sumabsu used: ", round(x$sumabsu,3), "\n")
  cat("Sumabsv used: ", round(x$sumabsv,3), "\n")
  if(x$upos) cat("U's constrained to be positive", fill=TRUE)
  if(x$uneg) cat("U's constrained to be negative", fill=TRUE)
  if(x$vpos) cat("V's constrained to be positive", fill=TRUE)
  if(verbose){
    for(k in 1:x$K){
      cat("\n Component ", k, ":\n")
      u <- x$u[,k]
      v <- x$v[,k]
      if(is.null(x$xnames)) x$xnames <- 1:length(u)
      if(is.null(x$znames)) x$znames <- 1:length(v)
      cat(fill=T)
      us <- cbind(x$xnames[u!=0], round(u[u!=0],3))
      dimnames(us) <- list(1:sum(u!=0), c("Row Feature Name", "Row Feature Weight"))
      vs <- cbind(x$znames[v!=0], round(v[v!=0],3))
      dimnames(vs) <- list(1:sum(v!=0), c("Column Feature Name", "Column Feature Weight"))
      print(us, quote=FALSE, sep="\t")
      cat(fill=T)
      print(vs, quote=FALSE, sep="\t")
    }  
  }
}


MultiSparseCCA <- function(x, z, K, sumabsu, sumabsv, v, niter=niter,trace=trace,upos,uneg,vpos){
  # do a bunch of single-factor sparse ccas; each time, do the next sparse cca on the residuals from the previous
  v.init <- v
  u=v=d=NULL
  xres <- x
  zres <- z
  for(k in 1:K){
    if(vpos && sum(abs(v.init[v.init[,k]>0,k]))<sum(abs(v.init[v.init[,k]<0,k]))) v.init[,k] <- -v.init[,k]
    out <- SparseCCA(xres,zres,sumabsu=sumabsu,sumabsv=sumabsv,niter=niter,v=v.init[,k],trace=trace,upos=upos,uneg=uneg,vpos=vpos)
    coef <- out$d 
    d <- c(d, coef)
    xres <- rbind(xres, sqrt(coef)*t(out$u))
    zres <- rbind(zres, -sqrt(coef)*t(out$v))
    u <- cbind(u, out$u)
    v <- cbind(v, out$v)
  }
  return(list(u=u,v=v,d=d,v.init=v.init, sumabsu=sumabsu, sumabsv=sumabsv, K=K,upos=upos,uneg=uneg,vpos=vpos))
}


