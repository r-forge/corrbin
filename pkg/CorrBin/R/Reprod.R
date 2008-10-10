
mc.est <- function(cbdata){
  #by trt
  do.est.fun <- function(x){
    est <- .Call("ReprodEstimates", as.integer(x$ClusterSize), as.integer(x$NResp), 
                   as.integer(x$Freq))
    est <- cbind(c(1,rep(NA,nrow(est)-1)), est) 
    idx <- upper.tri(est,diag=T)
    est.d <- data.frame(Prob=est[idx], ClusterSize=as.integer(col(est)[idx]-1), 
                        NResp=as.integer(row(est)[idx]-1),
                        Trt=x$Trt[1])
    est.d}  
  
  est.list <- by(cbdata, list(Trt=cbdata$Trt), do.est.fun)
  do.call(rbind,est.list)}

mc.test.chisq <- function(cbdata){
  cbdata <- subset(cbdata, Freq>0)
 
  get.T <- function(x){
      max.size <- max(x$ClusterSize)
      K <- sum(x$Freq)  
      K.r <- xtabs(Freq~ClusterSize, data=x)
      scores <- -(max.size - (2*(1:max.size)-1))/2
      a.r <- by(x, x$ClusterSize, function(z){sum(z$NResp*z$Freq)})[1:length(K.r)]
      cl.sizes <- as.numeric(names(a.r))
      N <- sum(cl.sizes*K.r)
      p.hat <- sum(a.r)/N
      p.hat.r <- a.r/(cl.sizes * K.r)
      T.stat <-  sum(cl.sizes * K.r * scores[cl.sizes]*p.hat.r)
      E.T.stat <- p.hat * sum(cl.sizes * K.r * scores[cl.sizes])
      rho.hat <- 1-sum(x$Freq*(x$ClusterSize-x$NResp)*x$NResp/x$ClusterSize)/((N-K)*p.hat*(1-p.hat))
      cat("rho=",rho.hat, "\n")
      Var0.T.stat <- p.hat*(1-p.hat)*sum(cl.sizes*K.r*scores[cl.sizes]^2*
                       (1+(cl.sizes-1)*rho.hat))
      b.r <- sum(cl.sizes * K.r * scores[cl.sizes])/N
      Var.E0T <- p.hat*(1-p.hat)*sum(cl.sizes*K.r*b.r^2*(1+(cl.sizes-1)*rho.hat))
      cov.T.E0T <- p.hat*(1-p.hat)*sum(cl.sizes*K.r*b.r*scores[cl.sizes]*
                  (1+(cl.sizes-1)*rho.hat))
      Var.T.stat <- Var0.T.stat + Var.E0T - 2*cov.T.E0T
      X.stat <- (T.stat - E.T.stat)^2/Var.T.stat
      X.stat}
      
   chis <- by(cbdata, cbdata$Trt, get.T)
   chis <- chis[1:length(chis)]
   chi.list <- list(chi.sq=chis, p=pchisq(chis, df=1, lower.tail=F))
   overall.chi <- sum(chis)
   overall.df <- length(chis)
   list(overall.chi=overall.chi, overall.p=pchisq(overall.chi, df=overall.df, lower.tail=F), 
        individual=chi.list)
}


mix.mc.mle <- function(cbdata, turn=1, control=mixControl()){ 
  attach(control)
  on.exit(detach(control))
  tab <- xtabs(Freq~factor(ClusterSize,levels=1:max(ClusterSize))+
                factor(NResp,levels=0:max(ClusterSize))+Trt, data=cbdata)
   size <- dim(tab)[1]
   ntrt <- dim(tab)[3]
   ntot <- sum(tab)
  storage.mode(tab) <- "double"
   Q <- array(0, dim=rep(size+1,ntrt))
   storage.mode(Q) <- "double"
     
   S <- DownUpMatrix(size, ntrt, turn)
   storage.mode(S) <- "integer"
   
   if ((start=="H0")&(method=="EM")){
     warning("The EM algorithm can only use 'start=uniform'. Switching options.")
     start <- "uniform"
  }
   if (start=="H0"){
     const.row <- matrix(0:size, nr=size+1, nc=ntrt)
     Q[const.row+1] <- 1/(size+1)
      }
   else {  #start=="uniform"
     Q[S+1] <- 1/(nrow(S))
    }
    
  res <- switch(method,
      EM = .Call("MixReprodQ", Q, S, tab, as.integer(max.iter), as.double(eps), 
                    as.integer(verbose)),
      ISDM = .Call("ReprodISDM", Q, S, tab, as.integer(max.iter), as.integer(max.directions),
                   as.double(eps),  as.integer(verbose)))
 
  names(res) <- c("MLest","Q","D","loglik", "converge")
 
  dimnames(res$MLest) <- list(NResp=0:size, ClusterSize=1:size, Trt=1:ntrt)
  res$MLest <- as.data.frame.table(res$MLest)
  names(res$MLest) <- c("NResp","ClusterSize","Trt","Prob") 
  res$MLest$NResp  <- as.numeric(as.character(res$MLest$NResp))
  res$MLest$ClusterSize  <- as.numeric(as.character(res$MLest$ClusterSize))
  res$MLest <- subset(res$MLest, NResp <= ClusterSize)
  levels(res$MLest$Trt) <- levels(cbdata$Trt)
  
  names(res$converge) <- c("rel.error", "n.iter")

  res
}


DownUpMatrix <- function(size, ntrt, turn){
  if ((turn<1)|(turn>ntrt)) stop("turn should be between 1 and ntrt")
  
    if (turn==1){
       res <- .Call("makeSmatrix", as.integer(size), as.integer(ntrt))
       return(res)
    }
    if (turn==ntrt){
       res <- .Call("makeSmatrix", as.integer(size), as.integer(ntrt))
       return(size - res)
    }
  
  
    res1 <- .Call("makeSmatrix", as.integer(size), as.integer(turn))
    res1 <- size - res1;
  
  res2list <- list()
  for (sq in 0:size){
    
      S <- .Call("makeSmatrix", as.integer(size-sq), as.integer(ntrt-turn))
      res2list <- c(res2list, list(sq+S))
    
  }
  
    res1list <- by(res1, res1[,turn], function(x)x)
    res <- mapply(merge, res1list, res2list, MoreArgs=list(by=NULL), SIMPLIFY=FALSE)
    res <- data.matrix(do.call(rbind, res))
    rownames(res) <- NULL
    colnames(res) <- NULL
    res
  
}

mixControl <- function(method=c("ISDM","EM"), eps=0.005, max.iter=5000, 
      max.directions=0, start=ifelse(method=="ISDM", "H0", "uniform"), verbose=FALSE){
  method <- match.arg(method)
  start <- match.arg(start, c("uniform","H0"))
  list(method = match.arg(method), eps = eps, max.iter = max.iter,
       max.directions = max.directions, start=start, verbose = verbose)
}

null.LRT <- function(cbdata, control=mixControl()){
   # LL under null hypothesis of equality (+ marginal compatibility)
   a <- with(cbdata, aggregate(Freq, list(ClusterSize=ClusterSize,NResp=NResp), sum))
   names(a)[names(a)=="x"] <- "Freq"
   a$ClusterSize <- as.integer(as.character(a$ClusterSize))
   a$NResp <- as.integer(as.character(a$NResp))
   a$Trt <- 1
                       
   b <- mc.est(a)
  b <- merge(cbdata, b, all.x=T, by=c("ClusterSize","NResp"))
  ll0 <- with(b, sum(Freq*log(Prob)))
  
  # LL under alternative hypothesis of stoch ordering (+ marginal compatibility)
  res <- mix.mc.mle(cbdata, control=control)
  ll1 <- res$loglik
  lrt <- 2*(ll1 - ll0)
  attr(lrt, "ll0") <- ll0
  attr(lrt, "ll1") <- ll1
  lrt
 }
  

perm.LRT <- function(cbdata, R=100, control=mixControl()){
   require(boot)
   dat2 <- cbdata[rep(1:nrow(cbdata), cbdata$Freq),]  #each row is one sample
   dat2$Freq <- NULL
   
   boot.LRT.fun <- function(dat, idx){
     dat.new <- cbind(dat[idx, c("ClusterSize","NResp")], Trt=dat$Trt)   #rearrange clusters
      dat.f <- aggregate(dat.new$Trt, 
                list(Trt=dat.new$Trt, ClusterSize=dat.new$ClusterSize, NResp=dat.new$NResp), length)
     names(dat.f)[names(dat.f)=="x"] <- "Freq"
    dat.f$ClusterSize <- as.numeric(as.character(dat.f$ClusterSize))
     dat.f$NResp <- as.numeric(as.character(dat.f$NResp))
                    
    stat <- null.LRT(dat.f, control=control)
    stat}        
        
   res <- boot(dat2, boot.LRT.fun, R=R, sim="permutation")
     
   p <- mean(res$t[,1] >= res$t0)
   LRT <- res$t0
   list(LRT=LRT, p.val=p, boot.res=res)}           

