
#'@import combinat

#'@rdname CorrBin-internal

    mChoose <- function(n, rvec, log=FALSE){
      rlast <- n - sum(rvec)
      rveclong <- c(rvec, rlast)
      if (any(rveclong < 0)) return(0)
      
      res <- lgamma(n + 1) - sum(lgamma(rveclong + 1))
      if (log) res else exp(res)
    }

tau <- function(cmdata, type=c("averaged","cluster")){
  type <- match.arg(type)
  
  
    nc <- attr(cmdata, "ncat")
    nrespvars <- paste("NResp", 1:nc, sep=".")
    M <- max(cmdata$ClusterSize)
  
  # multinomial lookup table
  mctab <- mChooseTable(M, nc, log=FALSE)
  
  res <- list()
  for (trt in levels(cmdata$Trt)){
    cm1 <- cmdata[cmdata$Trt==trt,]
    # observed freq lookup table
    atab <- array(0, dim=rep(M+1, nc))
    a.idx <- data.matrix(cm1[,nrespvars])
    atab[a.idx + 1] <- atab[a.idx + 1] + cm1$Freq
    
    if (type=="averaged"){
      Mn <- sum(cm1$Freq)
      
          
      res.trt <- array(NA, dim=rep(M+1, nc-1))
      dimnames(res.trt) <- rep.int(list(0:M), nc-1) 
      names(dimnames(res.trt)) <- paste("R", 1:(nc-1), sep="")
      # indices for possible values of r

         idx  <- hcube(rep( M +1,  nc-1 ))-1
          idxsum  <- rowSums(idx )
         idx  <- idx [ idxsum  <=  M , ,drop=FALSE]  #remove impossible indices
          idxsum  <-  idxsum [ idxsum  <=  M ]
      
      #indices for possible values of s 
      # (one more column than for r - ensures summation over all n's)

         sidx  <- hcube(rep( M +1,  nc ))-1
          sidxsum  <- rowSums(sidx )
         sidx  <- sidx [ sidxsum  <=  M , ,drop=FALSE]  #remove impossible indices
          sidxsum  <-  sidxsum [ sidxsum  <=  M ]
      
      for (i in 1:nrow(idx)){
        r <- idx[i,]
        s.idx <- which(sidxsum <= M-sum(r))
        lower.idx <- sidx[s.idx, , drop=FALSE]
        upper.idx <- lower.idx + rep(c(r,0), each=nrow(lower.idx))
        res.trt[rbind(r)+1] <- 
          sum(mctab[lower.idx+1] / mctab[upper.idx+1] * atab[upper.idx+1]) / Mn
      }
      
    } else {
      Mn <- xtabs(Freq ~ factor(ClusterSize, levels=1:M), data=cm1) 
      
      res.trt <- array(NA, dim=c(M, rep(M+1, nc-1))) #first dimension is 'n'
      dimnames(res.trt) <- c(list(1:M), rep.int(list(0:M), nc-1)) 
      names(dimnames(res.trt)) <- c("N",paste("R", 1:(nc-1), sep=""))
      for (n in which(Mn > 0)){
        # indices for possible values of r
        
           idx  <- hcube(rep( n +1,  nc-1 ))-1
            idxsum  <- rowSums(idx )
           idx  <- idx [ idxsum  <=  n , ,drop=FALSE]  #remove impossible indices
            idxsum  <-  idxsum [ idxsum  <=  n ]
        
        for (i in 1:nrow(idx)){
          r <- idx[i,]
          s.idx <- which(idxsum <= n-sum(r))
          lower.idx <- idx[s.idx, , drop=FALSE]
          upper.idx <- lower.idx + rep(r, each=nrow(lower.idx))
          lower.idx <- cbind(lower.idx, n-sum(r)-idxsum[s.idx])   #add implied last column
          upper.idx <- cbind(upper.idx, n-sum(r)-idxsum[s.idx])   #add implied last column
          res.trt[cbind(n,rbind(r)+1)] <- 
            sum(mctab[lower.idx+1] / mctab[upper.idx+1] * atab[upper.idx+1]) / Mn[n]
        }
      }
      
    }
    
    # append treatment-specific result to result list
    res.trt <- list(res.trt)
    names(res.trt) <- trt
    res <- c(res, res.trt) 
  }
  res
}

#'@rdname mc.est
#'@method mc.est CMData
#'@S3method mc.est CMData
#'@export
#'@param eps numeric; EM iterations proceed until the sum of squared changes fall below \code{eps}  


mc.est.CMData <- function(object, eps=1E-6, ...){

    nc <- attr(object, "ncat")      
    resp.vars1 <- paste("NResp", 1:(nc-1), sep=".")
   
    res <- mc.est.raw(object=object, eps=eps, ...)
    margres <- lapply(res, Marginals)  # has only NResp.1 - NResp.K
    
    mat.to.df <- function(idx, alist){
        dd <- as.data.frame.table(alist[[idx]], responseName="Prob")
        dd[c("N", resp.vars1)] <- lapply(dd[c("N", resp.vars1)], function(x)as.numeric(as.character(x)))
        dd$Trt <- names(alist)[idx]
        dd
    }
    margres <- lapply(1:length(margres), mat.to.df, alist=margres)
    fin <- do.call(rbind, margres)
    names(fin)[1] <- "ClusterSize"
    last.resp <- paste("NResp", nc, sep=".")
    fin[last.resp] <- fin$ClusterSize - rowSums(fin[resp.vars1]) # calculated omitted frequency
    fin$Trt <- factor(fin$Trt)
    fin[c("Trt","ClusterSize", resp.vars1, last.resp, "Prob")]
}

#'@rdname CorrBin-internal
Marginals <- function(theta){
  K <- length(dim(theta))
  M <- dim(theta)[1]-1
  
  res <- array(0, dim=c(M, rep(M+1, K)))
  dimnames(res) <- c(N=list(1:M), dimnames(theta))
  
  # indices for possible values of r
  
     idx  <- hcube(rep( M +1,  K+1 ))-1
      clustersize  <- rowSums(idx )
     idx  <- idx [ clustersize  <=  M , ,drop=FALSE]  #remove impossible indices
      clustersize  <-  clustersize [ clustersize  <=  M ]
  
  idx <- idx[ , -1, drop=FALSE]  #remove (K+1)st category
  
  
    curridx <- idx[clustersize==M, ,drop=FALSE]
    res[cbind(M, curridx+1)] <- theta[curridx+1]
  
  for (cs in seq.int(M-1,1)){
    
      curridx <- idx[clustersize==cs, , drop=FALSE]
      res[cbind(cs, curridx+1)] <- (cs+1- rowSums(curridx))/(cs+1) * res[cbind(cs+1, curridx+1)]
      for (j in 1:K){
        lookidx <- curridx
        lookidx[ ,j] <- lookidx[ ,j] + 1   #add 1 to the j-th coordinate
        res[cbind(cs, curridx+1)] <- res[cbind(cs, curridx+1)] + 
                                     lookidx[,j]/(cs+1) * res[cbind(cs+1, lookidx+1)]
      }  
    
  }
  
  res
}

#'@rdname CorrBin-internal
mc.est.raw <- function(object, ...) UseMethod("mc.est.raw")

#'@method mc.est.raw CMData
#'@S3method mc.est.raw CMData
mc.est.raw.CMData <- function(object, eps=1E-6, ...){
  cmdata <- object
  
    nc <- attr(cmdata, "ncat")
    nrespvars <- paste("NResp", 1:nc, sep=".")
    M <- max(cmdata$ClusterSize)
  
  
  # indices for possible values of r with clustersize = M
  
     idx  <- hcube(rep( M +1,  nc-1 ))-1
      idxsum  <- rowSums(idx )
     idx  <- idx [ idxsum  <=  M , ,drop=FALSE]  #remove impossible indices
      idxsum  <-  idxsum [ idxsum  <=  M ]
  

  res <- list()
  for (trt in levels(cmdata$Trt)){
    cm1 <- cmdata[cmdata$Trt==trt,]
    if (nrow(cm1) > 0){
      # observed freq lookup table
      atab <- array(0, dim=rep(M+1, nc))
      a.idx <- data.matrix(cm1[,nrespvars])
      atab[a.idx + 1] <- atab[a.idx + 1] + cm1$Freq
      Mn <- sum(cm1$Freq)
      
      
        res.trt <- array(NA, dim=rep(M+1, nc-1))
         
        #starting values
        res.trt[idx + 1] <- 1/nrow(idx)
        
        sqerror <- 1
        #EM update
        while (sqerror > eps){
              sqerror <- 0
              marg <- Marginals(res.trt)
          res.new <- array(NA, dim=rep(M+1, nc-1))
          res.new[idx + 1] <- 0
          
          
            for (i in 1:nrow(cm1)){
              rlong <- data.matrix(cm1[,nrespvars])[i,]    #nc elements
              r <- rlong[-nc]              #without the last category
              n <- cm1$ClusterSize[i]  
              # indices to which this cluster type contributes
              s.idx <- which(idxsum <= M-sum(r))
              tidx <- idx[s.idx, , drop=FALSE] + rep(r, each=length(s.idx))
              
              hvals <- apply(tidx, 1, function(tvec)prod(choose(tvec, r)) * choose(M-sum(tvec), n-sum(r))) 
              hvals <- hvals / choose(M, n)
              res.new[tidx+1] <- res.new[tidx+1] + atab[rbind(rlong)+1] / marg[rbind(c(n,r+1))] / Mn *
                                                   hvals * res.trt[tidx+1]
            }
          
              
          sqerror <- sum((res.new[idx+1] - res.trt[idx+1])^2)
              res.trt <- res.new 
        }
      
      
      # append treatment-specific result to result list
      dimnames(res.trt) <- rep.int(list(0:M), nc-1)
      names(dimnames(res.trt)) <- paste("NResp", 1:(nc-1), sep=".")
      res.trt <- list(res.trt)
    } else {
      res.trt <- list(c())
    } 
    res <- c(res, res.trt) 
  }
  names(res) <- levels(cmdata$Trt)
  res
} 
#'@rdname CorrBin-internal
tau.from.pi <- function(pimat){
  K <- length(dim(pimat))
  n <- dim(pimat)[1] - 1
  res <- array(NA, dim=rep(n+1, K)) 
  dimnames(res) <- rep.int(list(0:n), K) 
  names(dimnames(res)) <- paste("R", 1:K, sep="")

  # multinomial lookup table
  mctab <- mChooseTable(n, K+1, log=FALSE)
  
  # indices for possible values of r
  
     idx  <- hcube(rep( n +1,  K ))-1
      idxsum  <- rowSums(idx )
     idx  <- idx [ idxsum  <=  n , ,drop=FALSE]  #remove impossible indices
      idxsum  <-  idxsum [ idxsum  <=  n ]
  
  for (i in 1:nrow(idx)){
    r <- idx[i,]
    s.idx <- which(idxsum <= n-sum(r))
    lower.idx <- idx[s.idx, , drop=FALSE]
    upper.idx <- lower.idx + rep(r, each=nrow(lower.idx))
    lower.mc.idx <- cbind(lower.idx, n-sum(r)-idxsum[s.idx])   #add implied last column
    upper.mc.idx <- cbind(upper.idx, n-sum(r)-idxsum[s.idx])   #add implied last column
    res[rbind(r)+1] <- 
      sum(mctab[lower.mc.idx+1] / mctab[upper.mc.idx+1] * pimat[upper.idx+1])
  } 
  res
}

#'@rdname CorrBin-internal
p.from.tau <- function(taumat){
  K <- length(dim(taumat))
  idx <- diag(nrow=K)
  taumat[rbind(idx+1)]
}    

#'@rdname CorrBin-internal
corr.from.pi <- function(pimat){
  K <- length(dim(pimat))
  tt <- tau.from.pi(pimat)
  
  idx <- diag(nrow=K)
  numerator <- outer(1:K, 1:K, function(i,j){
     tt[idx[i,]+idx[j,]+1] - tt[idx[i,]+1] * tt[idx[j,]+1]})
  denominator <- outer(1:K, 1:K, function(i,j){
     tt[idx[i,]+1] * ifelse(i==j, 1-tt[idx[i,]+1], -tt[idx[j,]+1])})  
  res <- numerator / denominator    #the negative sign is in the denominator
}


#'@rdname mc.test.chisq
#'@method mc.test.chisq CMData
#'@S3method mc.test.chisq CMData
#'@export

mc.test.chisq.CMData <- function(object, ...){
  cmdata <- object[object$Freq > 0, ]
  K <- attr(object, "ncat")-1
  nrespvars <- paste("NResp", 1:K, sep=".")
  
  get.T <- function(x){
      x$Trt <- factor(x$Trt)  #remove unused levels
      pim <- mc.est.raw(x)[[1]]  #only one treatment group
      tt <- tau.from.pi(pim)
      p <- p.from.tau(tt)
      phi <- corr.from.pi(pim)
      xx <- x[rep(1:nrow(x), x$Freq),]
      xx$Freq <- 1
      
      M <- max(x$ClusterSize)
      Mn <- table(factor(xx$ClusterSize, levels=1:M)) 

      scores <- (1:M) - (M+1)/2
      
      Rmat <- data.matrix(xx[,nrespvars,drop=FALSE])
      nvec <- xx$ClusterSize
      cvec <- scores[nvec] 
      c.bar <- weighted.mean(cvec, w=nvec)
      cvec <- cvec - c.bar 
            
      X <- t(Rmat) %*% cvec
      Sigma <- diag(p, nrow=length(p)) - outer(p,p)  #multinomial vcov
      od.matrix <- matrix(0, nrow=K, ncol=K)  #over-dispersion matrix
      for (n in 1:M){
        od.matrix <- od.matrix + n * Mn[n] * (scores[n]-c.bar)^2 * (1+(n-1)*phi)
      }
      Sigma <- Sigma * od.matrix
      
      Tstat <- t(X) %*% solve(Sigma) %*% X       
      Tstat
   }
      
   chis <- by(cmdata, cmdata$Trt, get.T)
   chis <- chis[1:length(chis)]
   chi.list <- list(chi.sq=chis, p=pchisq(chis, df=K, lower.tail=FALSE))
   overall.chi <- sum(chis)
   overall.df <- length(chis) * K
   list(overall.chi=overall.chi, overall.p=pchisq(overall.chi, df=overall.df, lower.tail=FALSE), 
        individual=chi.list)
}

#'@rdname CorrBin-internal
  mChooseTable <- function(n, k, log=FALSE){
    res <- array(NA, dim=rep.int(n+1, k))
    dimnames(res) <- rep.int(list(0:n), k)
    
    idx <- hcube(rep.int(n+1, k)) - 1
    idx <- idx[rowSums(idx) <= n, ,drop=FALSE]
    for (i in 1:nrow(idx)){
        r <- idx[i, ]
        res[rbind(r)+1] <- mChoose(n=sum(r), rvec=r, log=log)
    }
    res
  }
