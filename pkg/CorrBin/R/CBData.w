\documentclass[reqno]{amsart}
\usepackage[margin=1in]{geometry}
\usepackage[colorlinks=true,linkcolor=blue]{hyperref}
\renewcommand{\NWtarget}[2]{\hypertarget{#1}{#2}}
\renewcommand{\NWlink}[2]{\hyperlink{#1}{#2}}
\newcommand{\bv}{\mathbf{v}}
\newcommand{\bq}{\mathbf{q}}
\newcommand{\bpi}{\text{\boldmath $\pi$}}
\newcommand{\leqst}{\mathrel{\preceq^{st}}}
\newcommand{\geqst}{\mathrel{\succeq^{st}}}

\title{Parametric analysis of correlated binary data}
\author{Aniko Szabo}
\date{\today}


\begin{document}
\begin{abstract} We define a class for describing data from toxicology experiments
and implement fitting of a variety of existing models and trend tests.
\end{abstract}
\maketitle

\section{Defining \texttt{CBData} -- a class for \textbf{C}lustered \textbf{B}inary \textbf{Data}}
We start with defining an S3 class describing data from toxicology experiments. The
class is a data frame with the following columns:

\begin{description}
\item[Trt] a factor defining (treatment) groups
\item[ClusterSize] an integer-valued variable defining the cluster size
\item[NResp] an integer-valued variable defining the number of responses (1s)
\item[Freq]  an integer-valued  variable defining frequency for each
Trt/ClusterSize/NResp combination
\end{description}

\texttt{CBData} converts a data frame to a CBData object. \texttt{x}
is the input data frame; \texttt{trt}, \texttt{clustersize}, \texttt{nresp} and
\texttt{freq} could be strings or column indices defining the appropriate
variable in \texttt{x} (\texttt{freq} can also be NULL, in which case it is 
assumed that each combination has frequency 1).
@O CBData.R @{
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
@| CBdata.data.frame  @}


The \texttt{read.CBData} function reads in clustered binary data from a tab-delimited
text file. The first column should give the treatment group, the second the size of the cluster,
the third the number of responses in the cluster. Optionally, a fourth column could
give the number of times the given combination occurs in the data.

@o CBData.R
@{
read.CBData <- function(file, with.freq=TRUE, ...){
  d <- read.table(file, col.names=c("Trt","ClusterSize","NResp", if (with.freq) "Freq"), ...)
  if (!with.freq) d$Freq <- 1
  d <- aggregate(d$Freq, list(Trt=d$Trt, ClusterSize=d$ClusterSize, NResp=d$NResp),sum)
  names(d)[4] <- "Freq"
  d$ClusterSize <- as.numeric(as.character(d$ClusterSize))
  d$NResp <- as.numeric(as.character(d$NResp))
  d <- CBData(d, "Trt", "ClusterSize", "NResp", "Freq")
  d}
@| read.CBdata @}

\texttt{unwrap.CBData} is a utility function that reformats a CBData object so that
each row is one observation (instead of one cluster). A new `ID' variable is added
to indicate clusters. 

@O CBData.R @{

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
@| unwrap.CBData @}

\section{Rao-Scott adjusted Cochran-Armitage test}
The RS-adjusted CA test for trend is based on design-effect adjustment.

@O CBData.R @{
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
  p.val <- pnorm(RS, lower.tail=F)
  list(statistic=RS, p.val=p.val)
  }
@| RS.trend.test @}

\section{GEE based test}

@O CBData.R @{
require(geepack)

GEE.trend.test <- function(cbdata, scale.method=c("fixed", "trend", "all")){
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
  p <- pnorm(testres, lower.tail=F)
  list(statistic=testres, p.val=p)
 }   
@| GEE.trend.test @} 

\section{Generating random data}
\texttt{ran.mc.CBData} generates a random CBData object from a given two-parameter
distribution. \texttt{sample.sizes} is a dataset with variables Trt, ClusterSize and
Freq giving the number of clusters to be generated for each Trt/ClusterSize combination.
\texttt{p.gen.fun} and \texttt{rho.gen.fun} are functions that generate the parameter
values for each treatment group ($g=1$ corresponds to the lowest group, $g=2$ to the
second, etc). \texttt{pdf.fun} is a function(p, rho, n) generating the pdf of the
number of responses given the two parameters \texttt{p} and \texttt{rho}, and the
cluster size \texttt{n}.

@O CBData.R @{
ran.mc.CBData <- function(sample.sizes, p.gen.fun=function(g)0.3,
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

   a <- by(sample.sizes, list(Trt=sample.sizes$Trt, ClusterSize=sample.sizes$ClusterSize), ran.gen)
   a <- data.frame(do.call(rbind, a))
   a$Trt <- factor(a$Trt, labels=levels(sample.sizes$Trt))
   a <- subset(a, Freq>0)
   class(a) <-  c("CBData", "data.frame")
   a
 }
@| ran.mc.CBData @}

\subsection{Parametric pdf generating functions}
\texttt{betabin.pdf} and \texttt{qpower.pdf} provide two classic distributions --
beta-binomial and q-power -- for generating correlated binary data. Either
can be used in \texttt{ran.mc.CBData}.
@o CBData.R
@{
 betabin.pdf <- function(p, rho, n){
   a <- p*(1/rho-1)
   b <- (1-p)*(1/rho-1)
   idx <- 0:n
   res <- choose(n, idx)*beta(a+idx, b+n-idx)/beta(a,b)
   res
  } 
@}
@o CBData.R
@{
 qpower.pdf <- function(p, rho, n){
   .q <- 1-p
   gamm <- log2(log(.q^2+rho*.q*(1-.q))/log(.q))
   res <- numeric(n+1)
   for (y in 0:n){
     idx <- 0:y
     res[y+1] <- choose(n,y) * sum( (-1)^idx * choose(y,idx) * .q^((n-y+idx)^gamm))
   }
   #if (any(res<0, res>1)) print(c(p,rho))
   res <- pmax(pmin(res,1),0)  #to account for numerical imprecision
   res
 }
@}
\section{Fitting parametric dose-response models}
To provide flexibility with selecting the parametric family we follow the approach 
of \texttt{glm} in defining the appropriate lists containing the log-likelihood and
its gradient.

\subsection{Beta-binomial model}

The beta-binomial distribution is parameterized by the marginal mean $p$ and 
the intra-cluster correlation $r$. \texttt{X.p} and \texttt{X.r} are the 
design matrices for modeling them as a function of the group, \texttt{param.p} and 
\texttt{param.r} are the corresponding parameters.

@o CBData.R @{
betabinomial <- function(link.p="logit", link.r="logit"){ 
  fam.name <- "beta-binomial"
  @< Set up p-link function @>
  @< Set up r-link function @>
  @< Define betabinomial log-likelihood @>
  @< Define betabinomial gradient @>


  res <- list(name=fam.name, logl = logl, gr=gr)  
  class(res) <- "CBfamily"
  res}
@| betabinomial @}
  
The code for setting up the link function is taken from the \texttt{binomial} 
function.
@D Set up p-link function @{
    linktemp.p <- substitute(link.p)
    if (!is.character(linktemp.p)) {
        linktemp.p <- deparse(linktemp.p)
    }
    okLinks <- c("logit", "probit", "cloglog", "cauchit")
    if (linktemp %in% okLinks) 
        stats <- make.link(linktemp.p)
    else if (is.character(link.p)) {
        stats <- make.link(link.p)
        linktemp.p <- link.p
    }
    else {
        if (inherits(link.p, "link-glm")) {
            stats <- link.p
            if (!is.null(stats$name)) 
                linktemp.p <- stats$name
        }
        else {
            stop(gettextf("link \"%s\" not available for the %s family; available links are %s", 
                linktemp.p, fam.name, paste(sQuote(okLinks), collapse = ", ")), 
                domain = NA)
        }
    }
 
@}

For the correlation coefficient the link function is actually for $2*(r+1)$. This way
the usual links for 0 -- 1 data apply to the $-1$ -- 1 range of $r$.

@D Set up r-link function @{
   linktemp.r <- substitute(link.r)
    if (!is.character(linktemp.r)) {
        linktemp.r <- deparse(linktemp.r)
    }
    okLinks <- c("logit", "probit", "cloglog", "cauchit")
    if (linktemp %in% okLinks) 
        stats <- make.link(linktemp.r)
    else if (is.character(link.r)) {
        stats <- make.link(link.r)
        linktemp.r <- link.r
    }
    else {
        if (inherits(link.r, "link-glm")) {
            stats <- link.r
            if (!is.null(stats$name)) 
                linktemp.r <- stats$name
        }
        else {
            stop(gettextf("link \"%s\" not available for the %s family; available links are %s", 
                linktemp.r, fam.name, paste(sQuote(okLinks), collapse = ", ")), 
                domain = NA)
        }
    }
@}

The beta-binomial log-likelihood is
\begin{equation}
\log L_{bb}(x)=
\end{equation}

@D Define betabinomial log-likelihood @{
  logl <- function(param.p, param.r, X.p, X.r){

  b.p<-param[param.ind=="p"]
  b.r<-param[param.ind=="r"]
  
  score.p <- Z%*%b.p
  score.r <- Z%*%b.r
  
  p <- inv.logit(score.p)
  r <- inv.logit(score.r)
  
  # compute log-Likelihood
  logL <- 0  
  for (i in 1:length(ClusterSize))
    {
      n <- ClusterSize[i]
      y <- NResp[i]
      vec1<-0:(y-1)
      vec2<-0:(n-y-1)
      vec3<-0:(n-1)
      l<-ifelse(y>0,sum(log(p[i]+r[i]*vec1)),0) + ifelse(y<n,sum(log(1-p[i]+r[i]*vec2)),0) - sum(log(1+r[i]*vec3)) 
    logL <- logL+l    
    }
  return(-logL)
}

@}

\end{document}