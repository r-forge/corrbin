
CBData <- function(x, trt, clustersize, nresp, freq=NULL){
  if (!is.data.frame(x)) stop("x has to be a data frame")
  nms <- names(x)
  process.var <- function(var){
    if (is.character(var)){
       if (var %in% nms) res <- x[[var]]
       else stop(paste("Variable '", var, "' not found"))
    }
    else {
      if (is.numeric(var)){
         if (var %in% seq(along=nms)) res <- x[[var]]
         else stop(paste("Column", var, " not found"))
      }
      else stop(paste("Invalid variable specification:",var))
    }
  }
  trtvar <- factor(process.var(trt))
  csvar <- process.var(clustersize)
  nrespvar <- process.var(nresp)
  if (is.null(freq)) freqvar <- rep(1, nrow(x))
  else freqvar <- process.var(freq)
  
  d <- data.frame(Trt=trtvar, ClusterSize=csvar, NResp=nrespvar, Freq=freqvar)
  d <- aggregate(d$Freq, list(Trt=d$Trt, ClusterSize=d$ClusterSize, NResp=d$NResp),sum)
  names(d)[4] <- "Freq"
  d$ClusterSize <- as.numeric(as.character(d$ClusterSize))
  d$NResp <- as.numeric(as.character(d$NResp))
  class(d) <- c("CBData", "data.frame")
  d}

read.CBData <- function(file, with.freq=TRUE, ...){
  d <- read.table(file, col.names=c("Trt","ClusterSize","NResp", if (with.freq) "Freq"), ...)
  if (!with.freq) d$Freq <- 1
  d <- aggregate(d$Freq, list(Trt=d$Trt, ClusterSize=d$ClusterSize, NResp=d$NResp),sum)
  names(d)[4] <- "Freq"
  d$ClusterSize <- as.numeric(as.character(d$ClusterSize))
  d$NResp <- as.numeric(as.character(d$NResp))
  d <- CBData(d, "Trt", "ClusterSize", "NResp", "Freq")
  d}


unwrap.CBData <- function(cbdata){
  freqs <- rep(1:nrow(cbdata), cbdata$Freq)
  cb1 <- cbdata[freqs,]
  cb1$Freq <- NULL
  cb1$ID <- factor(1:nrow(cb1))
  pos.idx <- rep(1:nrow(cb1), cb1$NResp)
  cb.pos <- cb1[pos.idx,]
  cb.pos$Resp <- 1
  cb.pos$NResp <- NULL
  neg.idx <- rep(1:nrow(cb1), cb1$ClusterSize-cb1$NResp)
  cb.neg <- cb1[neg.idx,]
  cb.neg$Resp <- 0
  cb.neg$NResp <- NULL
  res <- rbind(cb.pos, cb.neg)
  res[order(res$ID),]
  }

RS.trend.test <- function(cbdata){  
   dat2 <- cbdata[rep(1:nrow(cbdata), cbdata$Freq),]  #each row is one sample
   dat2$Trt <- factor(dat2$Trt)  #remove unused levels
  attach(dat2)
  on.exit(detach(dat2))
  x.i <- pmax(tapply(NResp, Trt, sum), 0.5)  #"continuity" adjustment to avoid RS=NaN
  n.i <- tapply(ClusterSize, Trt, sum)
  m.i <- table(Trt)
  p.i.hat <- x.i/n.i
  r.ij <- NResp - ClusterSize*p.i.hat[Trt]
  v.i <- m.i/(m.i-1)/n.i^2*tapply(r.ij^2, Trt, sum)
  d.i <- n.i * v.i / (p.i.hat*(1-p.i.hat))   #design effect
  x.i.new <- x.i/d.i
  n.i.new <- n.i/d.i
  p.hat <- sum(x.i.new)/sum(n.i.new)
  
  scores <- (1:nlevels(Trt))-1
  mean.score <- sum(scores*n.i.new)/sum(n.i.new)
  var.scores <- sum(n.i.new*(scores-mean.score)^2)
  RS <- (sum(x.i.new*scores) - p.hat*sum(n.i.new*scores)) / 
        sqrt(p.hat*(1-p.hat)*var.scores)
  p.val <- pnorm(RS, lower.tail=FALSE)
  list(statistic=RS, p.val=p.val)
  }


GEE.trend.test <- function(cbdata, scale.method=c("fixed", "trend", "all")){
  require(geepack)
  ucb <- unwrap.CBData(cbdata)
  scale.method <- match.arg(scale.method)
  if (scale.method=="fixed") {
    geemod <- geese(Resp~unclass(Trt), id=ID, scale.fix=FALSE, data=ucb,
                    family=binomial, corstr="exch") }  
  else if (scale.method=="trend"){
    geemod <- geese(Resp~unclass(Trt), sformula=~unclass(Trt), id=ID,  data=ucb,
                   family=binomial, sca.link="log", corstr="exch")}
  else if (scale.method=="all"){
    geemod <- geese(Resp~unclass(Trt), id=ID,  sformula=~Trt, data=ucb,
                    family=binomial, sca.link="log", corstr="exch") } 
  geesum <- summary(geemod)
  testres <- geesum$mean[2,"estimate"]/geesum$mean[2,"san.se"]
  p <- pnorm(testres, lower.tail=FALSE)
  list(statistic=testres, p.val=p)
 }   

ran.CBData <- function(sample.sizes, p.gen.fun=function(g)0.3,
                           rho.gen.fun=function(g)0.2, pdf.fun=qpower.pdf){
   ran.gen <- function(d){
   # d is subset(sample.sizes, Trt==trt, ClusterSize==cs)
     cs <- d$ClusterSize[1]
     trt <- unclass(d$Trt)[1]
     n <- d$Freq[1]
     p <- p.gen.fun(trt)
     rho <- rho.gen.fun(trt)
     probs <- pdf.fun(p, rho, cs)
     tmp <- rmultinom(n=1, size=n, prob=probs)[,1]
     cbind(Freq=tmp, NResp=0:cs, ClusterSize=d$ClusterSize, Trt=d$Trt)}

   sst <- if (is.factor(sample.sizes$Trt)) sample.sizes$Trt else factor(sample.sizes$Trt)
   a <- by(sample.sizes, list(Trt=sst, ClusterSize=sample.sizes$ClusterSize), ran.gen)
   a <- data.frame(do.call(rbind, a))
   a$Trt <- factor(a$Trt, labels=levels(sst))
   a <- subset(a, Freq>0)
   class(a) <-  c("CBData", "data.frame")
   a
 }

 betabin.pdf <- function(p, rho, n){
   a <- p*(1/rho-1)
   b <- (1-p)*(1/rho-1)
   idx <- 0:n
   res <- choose(n, idx)*beta(a+idx, b+n-idx)/beta(a,b)
   res
  } 

 qpower.pdf <- function(p, rho, n){
   .q <- 1-p
   gamm <- log2(log(.q^2+rho*.q*(1-.q))/log(.q))
   res <- numeric(n+1)
   for (y in 0:n){
     idx <- 0:y
     res[y+1] <- choose(n,y) * sum( (-1)^idx * choose(y,idx) * .q^((n-y+idx)^gamm))
   }
   res <- pmax(pmin(res,1),0)  #to account for numerical imprecision
   res
 }
