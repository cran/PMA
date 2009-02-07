CGH.SparseCCA <- function(x,y,lam1,lam2,sumabsu,chrom=NULL,niter=20,v=NULL,trace=trace, upos, uneg){
  # x is gene expression data, y is cgh data
  vold <- rnorm(ncol(y))
  for(i in 1:niter){
    if(sum(abs(vold-v))>1e-5 && mean(v==0)<1){
      vold <- v
      if(trace) cat(i,fill=FALSE)
      argu <- t(x)%*%(y%*%v)
      if(upos) argu <- pmax(0, argu)
      if(uneg) argu <- pmin(0, argu)
      lamu <- BinarySearch(argu,sumabsu)
      u <- soft(argu,lamu)/l2n(soft(argu,lamu))
      vnew <- numeric(ncol(y))
      xu <- x%*%u
      for(j in sort(unique(chrom))){
        yxu <- as.numeric(t(y[,chrom==j])%*%xu)
        coefs <- FLSA(yxu/l2n(yxu),lambda1=lam1,lambda2=lam2)
        vnew[chrom==j] <- coefs
      }    
      v <- vnew
      if(sum(abs(v))>0) v <- v/l2n(v)
      v <- matrix(v,ncol=1)
    }
  }
  if(mean(v==0)==1) u <- matrix(0, nrow=ncol(x), ncol=1)
  return(list(u=u,v=v))
}

CCAL1FL <- function(x,z,K=1,sumabsu=5,lambda=NULL,chrom=NULL,niter=20,v=NULL, trace=TRUE, standardize=TRUE, xnames=NULL, znames=NULL, upos=FALSE, uneg=FALSE){
  call <- match.call()
  if(nrow(x)!=nrow(z)) stop("x and z must have same number of rows")
  if(standardize){
    x <- scale(x, TRUE, TRUE)
    z <- scale(z, TRUE, TRUE)
  }
  if(sumabsu<1 || sumabsu>sqrt(ncol(x))) stop("sumabsu must be between 1 and sqrt(ncol(x))")
  v <- CheckVs(v,x,z,K)
  out <- MultiCGH.SparseCCA(x=x,y=z,K=K,sumabsu=sumabsu,lambda=lambda,chrom=chrom,niter=niter, trace=trace,v=v, upos=upos, uneg=uneg)
  returnobj <-  list(u=out$u,v=out$v,v.init=out$v.init, d=out$d, sumabsu=sumabsu, lambda=out$lambda, call=call, xnames=xnames, znames=znames, K=K, upos=upos, uneg=uneg)
  class(returnobj) <- "ccal1fl"
  return(returnobj)
}

print.ccal1fl <- function(x,verbose=FALSE,...){
  cat("Call: ")
  dput(x$call)
  cat("\n\n")
  cat("Num non-zero u's: ", apply(x$u!=0,2,sum),"\n")
  cat("Num non-zero v's: ", apply(x$v!=0,2,sum),"\n")
  cat("Sumabsu: ", x$sumabsu, "\n")
  cat("Lambda: ", x$lambda, "\n")
  if(x$upos) cat("U's constrained to be positive.", fill=TRUE)
  if(x$uneg) cat("U's constrained to be negative.", fill=TRUE)
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

MultiCGH.SparseCCA <- function(x,y,K,sumabsu,lambda,chrom,niter,v, trace=trace,upos,uneg){
  if(is.null(chrom)) chrom <- rep(1,ncol(y))
  uans = vans = dans = NULL
  if(is.null(lambda)) lambda <- ChooseLambda1Lambda2(as.numeric(v[,1]))
  xres <- x
  yres <- y
  uans <- matrix(NA, nrow=ncol(x), ncol=K)
  vans <- matrix(NA, nrow=ncol(y), ncol=K)
  dans <- numeric(K)
  for(k in 1:K){
    if(trace) cat(k,fill=TRUE)
    out <- CGH.SparseCCA(xres,yres,lambda,lambda,sumabsu,chrom=chrom,niter=niter,v=v[,k],trace=trace,upos=upos,uneg=uneg)
    if(sum(is.na(out$u))>0 || sum(is.na(out$v))>0 || mean(out$u==0)==1 || mean(out$v==0)==1){
      out$u <- matrix(0,nrow=ncol(x),ncol=1)
      out$v <- matrix(0,nrow=ncol(y),ncol=1)
    }
    uans[,k] <- out$u
    vans[,k] <- out$v
    coef <- sum((x%*%out$u)*(y%*%out$v))
    dans[k] <- coef
    xres <- rbind(xres, sqrt(coef)*t(out$u))
    yres <- rbind(yres, -sqrt(coef)*t(out$v))
  }
  if(trace) cat(fill=TRUE)
  return(list(u=uans,v=vans, dans=dans, lambda=lambda,v.init=v))
}


