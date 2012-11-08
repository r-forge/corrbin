
mc.est <- function(cbdata){
  #by trt
  do.est.fun <- function(x){
    est <- .Call("ReprodEstimates", as.integer(x$ClusterSize), as.integer(x$NResp), 
                   as.integer(x$Freq),PACKAGE="CorrBin")
    est <- cbind(c(1,rep(NA,nrow(est)-1)), est) 
    idx <- upper.tri(est,diag=TRUE)
    est.d <- data.frame(Prob=est[idx], ClusterSize=as.integer(col(est)[idx]-1), 
                        NResp=as.integer(row(est)[idx]-1),
                        Trt=x$Trt[1])
    est.d}  
  
  est.list <- by(cbdata, list(Trt=cbdata$Trt), do.est.fun)
  do.call(rbind,est.list)}

  
mc.test.chisq <- function(cbdata){
  cbdata <- cbdata[cbdata$Freq>0, ]
 
  get.T <- function(x){
      max.size <- max(x$ClusterSize)
      scores <- (1:max.size) - (max.size+1)/2
      p.hat <- with(x, sum(Freq*NResp) / sum(Freq*ClusterSize))
      rho.hat <- with(x, 1-sum(Freq*(ClusterSize-NResp)*NResp/ClusterSize) / 
          (sum(Freq*(ClusterSize-1))*p.hat*(1-p.hat)))  #Fleiss-Cuzick estimate
      c.bar <- with(x, sum(Freq*scores[ClusterSize]*ClusterSize) / sum(Freq*ClusterSize))
      T.center <- with(x, sum(Freq*(scores[ClusterSize]-c.bar)*NResp))
      Var.T.stat <-  with(x, 
         p.hat*(1-p.hat)*sum(Freq*(scores[ClusterSize]-c.bar)^2*ClusterSize*(1+(ClusterSize-1)*rho.hat)))
      X.stat <- (T.center)^2/Var.T.stat
      X.stat}
      
   chis <- by(cbdata, cbdata$Trt, get.T)
   chis <- chis[1:length(chis)]
   chi.list <- list(chi.sq=chis, p=pchisq(chis, df=1, lower.tail=FALSE))
   overall.chi <- sum(chis)
   overall.df <- length(chis)
   list(overall.chi=overall.chi, overall.p=pchisq(overall.chi, df=overall.df, lower.tail=FALSE), 
        individual=chi.list)
}


SO.mc.est <- function(cbdata, turn=1, control=soControl()){ 
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
   
   if ((control$start=="H0")&(control$method=="EM")){
     warning("The EM algorithm can only use 'start=uniform'. Switching options.")
     start <- "uniform"
  }
   if (control$start=="H0"){
     const.row <- matrix(0:size, nrow=size+1, ncol=ntrt)
     Q[const.row+1] <- 1/(size+1)
      }
   else {  #start=="uniform"
     Q[S+1] <- 1/(nrow(S))
    }
    
  res0 <- switch(control$method,
      EM = .Call("MixReprodQ", Q, S, tab, as.integer(control$max.iter), as.double(control$eps), 
                    as.integer(control$verbose), PACKAGE="CorrBin"),
      ISDM = .Call("ReprodISDM", Q, S, tab, as.integer(control$max.iter), as.integer(control$max.directions),
                   as.double(control$eps),  as.integer(control$verbose), PACKAGE="CorrBin"))
 
  names(res0) <- c("MLest","Q","D","loglik", "converge")
  names(res0$converge) <- c("rel.error", "n.iter")
  res <- res0$MLest
 
  dimnames(res) <- list(NResp=0:size, ClusterSize=1:size, Trt=1:ntrt)
  res <- as.data.frame.table(res)
  names(res) <- c("NResp","ClusterSize","Trt","Prob") 
  res$NResp  <- as.numeric(as.character(res$NResp))
  res$ClusterSize  <- as.numeric(as.character(res$ClusterSize))
  res <- res[res$NResp <= res$ClusterSize,]
  levels(res$Trt) <- levels(cbdata$Trt)
  
  attr(res, "loglik") <- res0$loglik
  attr(res, "converge") <- res0$converge
  res
}

soControl <- function(method=c("ISDM","EM"), eps=0.005, max.iter=5000, 
      max.directions=0, start=ifelse(method=="ISDM", "H0", "uniform"), verbose=FALSE){
  method <- match.arg(method)
  start <- match.arg(start, c("uniform","H0"))
  list(method = match.arg(method), eps = eps, max.iter = max.iter,
       max.directions = max.directions, start=start, verbose = verbose)
}

DownUpMatrix <- function(size, ntrt, turn){
  if ((turn<1)|(turn>ntrt)) stop("turn should be between 1 and ntrt")
  
    if (turn==1){
       res <- .Call("makeSmatrix", as.integer(size), as.integer(ntrt),PACKAGE="CorrBin")
       return(res)
    }
    if (turn==ntrt){
       res <- .Call("makeSmatrix", as.integer(size), as.integer(ntrt),PACKAGE="CorrBin")
       return(size - res)
    }
  
  
    res1 <- .Call("makeSmatrix", as.integer(size), as.integer(turn),PACKAGE="CorrBin")
    res1 <- size - res1;
  
  res2list <- list()
  for (sq in 0:size){
    
      S <- .Call("makeSmatrix", as.integer(size-sq), as.integer(ntrt-turn),PACKAGE="CorrBin")
      res2list <- c(res2list, list(sq+S))
    
  }
  
    res1list <- by(res1, res1[,turn], function(x)x)
    res <- mapply(merge, res1list, res2list, MoreArgs=list(by=NULL), SIMPLIFY=FALSE)
    res <- data.matrix(do.call(rbind, res))
    rownames(res) <- NULL
    colnames(res) <- NULL
    res
  
}

SO.LRT <- function(cbdata, control=soControl()){
   # LL under null hypothesis of equality (+ reproducibility)
   a <- with(cbdata, aggregate(Freq, list(ClusterSize=ClusterSize,NResp=NResp), sum))
   names(a)[names(a)=="x"] <- "Freq"
   a$ClusterSize <- as.integer(as.character(a$ClusterSize))
   a$NResp <- as.integer(as.character(a$NResp))
   a$Trt <- 1
                       
   b <- mc.est(a)
  b <- merge(cbdata, b, all.x=TRUE, by=c("ClusterSize","NResp"))
  ll0 <- with(b, sum(Freq*log(Prob)))
  
  # LL under alternative hypothesis of stoch ordering (+ reproducibility)
  res <- SO.mc.est(cbdata, control=control)
  ll1 <- attr(res, "loglik")
  lrt <- 2*(ll1 - ll0)
  attr(lrt, "ll0") <- ll0
  attr(lrt, "ll1") <- ll1
  lrt
 }
  

SO.trend.test <- function(cbdata, R=100, control=soControl()){
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
                    
    stat <- SO.LRT(dat.f, control=control)
    stat}        
        
   res <- boot(dat2, boot.LRT.fun, R=R, sim="permutation")
     
   p <- mean(res$t[,1] >= res$t0)
   LRT <- res$t0
   list(LRT=LRT, p.val=p, boot.res=res)}           

trend.test <- function(cbdata, test=c("RS","GEE","GEEtrend","GEEall","SO"), exact=test=="SO", 
                       R=100, control=soControl()){ 
   test <- match.arg(test)
   if (!exact && !(test=="SO")){
     res <- switch(test, RS=RS.trend.test(cbdata), 
                         GEE=GEE.trend.test(cbdata,scale.method="fixed"),
                         GEEtrend=GEE.trend.test(cbdata,scale.method="trend"),
                         GEEall=GEE.trend.test(cbdata,scale.method="all"))
   }
   else {
     dat2 <- cbdata[rep(1:nrow(cbdata), cbdata$Freq),]  #each row is one sample
     dat2$Freq <- NULL
     
     boot.LRT.fun <- function(dat, idx){
       dat.new <- cbind(dat[idx, c("ClusterSize","NResp")], Trt=dat$Trt)   #rearrange clusters
       dat.f <- aggregate(dat.new$Trt, 
                  list(Trt=dat.new$Trt, ClusterSize=dat.new$ClusterSize, NResp=dat.new$NResp), length)
       names(dat.f)[names(dat.f)=="x"] <- "Freq"
       dat.f$ClusterSize <- as.numeric(as.character(dat.f$ClusterSize))
       dat.f$NResp <- as.numeric(as.character(dat.f$NResp))
                      
       stat <- switch(test, SO=SO.LRT(dat.f, control=control),
                            RS=RS.trend.test(dat.f)$statistic,
                            GEE=GEE.trend.test(dat.f, scale.method="fixed")$statistic,
                            GEEtrend=GEE.trend.test(cbdata,scale.method="trend")$statistic,
                            GEEall=GEE.trend.test(cbdata,scale.method="all")$statistic)
       stat}        
          
     bootres <- boot(dat2, boot.LRT.fun, R=R, sim="permutation")
     res <- list(statistic=bootres$t0, p.val= mean(bootres$t[,1] >= bootres$t0))
     attr(res, "boot") <- bootres
   }   
   res}

NOSTASOT <- function(cbdata, test=c("RS","GEE","GEEtrend","GEEall","SO"), exact=test=="SO",
                     R=100, sig.level=0.05, control=soControl()){
   ntrt <- nlevels(cbdata$Trt)
   control.gr <- levels(cbdata$Trt)[1]
   p.vec <- array(NA, ntrt-1)
   names(p.vec) <- levels(cbdata$Trt)[-1]
   NOSTASOT.found <- FALSE
   curr.gr.idx <- ntrt
   curr.gr <- levels(cbdata$Trt)[ntrt]
   
   while (!NOSTASOT.found & (curr.gr.idx>1)){
     d1 <- cbdata[unclass(cbdata$Trt)<=curr.gr.idx, ]
     d1$Trt <- factor(d1$Trt) #eliminate unused levels
     tr.test <- trend.test(d1, test=test, exact=exact, R=R, control=control)
     p.vec[curr.gr] <- tr.test$p.val
     if (tr.test$p.val < sig.level){ #NOSTASOT not found yet
       curr.gr.idx <- curr.gr.idx - 1
       curr.gr <- levels(cbdata$Trt)[curr.gr.idx]
     }
     else { #NOSTASOT
       NOSTASOT.found <- TRUE
     }
   }
       
   list(NOSTASOT = curr.gr, p=p.vec)    
}       
