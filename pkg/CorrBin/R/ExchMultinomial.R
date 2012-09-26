
  library(combinat)


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
    cm1 <- subset(cmdata, Trt==trt)
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

mc.est <- function(cmdata, eps=1E-6){
  
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
    cm1 <- subset(cmdata, Trt==trt)
    # observed freq lookup table
    atab <- array(0, dim=rep(M+1, nc))
    a.idx <- data.matrix(cm1[,nrespvars])
    atab[a.idx + 1] <- atab[a.idx + 1] + cm1$Freq
    Mn <- sum(cm1$Freq)
    
    
      res.trt <- array(NA, dim=rep(M+1, nc-1))
      dimnames(res.trt) <- rep.int(list(0:M), nc-1)
       
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
    res.trt <- list(res.trt)
    names(res.trt) <- trt
    res <- c(res, res.trt) 
  }
  res
}
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
