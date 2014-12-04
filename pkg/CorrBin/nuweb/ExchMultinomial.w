\documentclass[reqno]{amsart}
\usepackage[margin=1in]{geometry}
\usepackage[colorlinks=true,linkcolor=blue]{hyperref}
\renewcommand{\NWtarget}[2]{\hypertarget{#1}{#2}}
\renewcommand{\NWlink}[2]{\hyperlink{#1}{#2}}

\providecommand{\tsum}{\textstyle\sum}
\providecommand{\dsum}{\displaystyle\sum}
\newcommand{\rvec}{\mathbf{r}}
\newcommand{\svec}{\mathbf{s}}
\newcommand{\tvec}{\mathbf{t}}
\newcommand{\dvec}{\mathbf{d}}
\newcommand{\taur}[1]{\tau_{r_1,\ldots,r_{#1}}}
\newcommand{\taurn}[2]{\tau_{r_1,\ldots,r_{#1}|#2}}
\newcommand{\htaurn}[2]{\hat{\tau}_{r_1,\ldots,r_{#1}|#2}}
\newcommand{\thetar}[1]{\theta_{r_1,\ldots,r_{#1}}}
\newcommand{\hthetar}[1]{\hat{\theta}_{r_1,\ldots,r_{#1}}}
\newcommand{\Arn}[2]{A_{r_1,\ldots,r_{#1}|#2}}
\newcommand{\X}{\mathcal{X}}
\newcommand{\V}{\mathcal{V}}

\DeclareMathOperator{\Prob}{P}

\title{Exchangeable model for multinomial data}
\author{Aniko Szabo}
\date{\today}

\begin{document}
\maketitle


\begin{abstract}
We implement parameter estimation for exchangeable multinomial data, including estimation under marginal compatibility. 
\end{abstract}

\section{Preliminaries}

We will be using object of \texttt{CMData} class, which is defined in \texttt{CMData.w}.

We will also need to load support libraries.

@o ..\R\ExchMultinomial.R @{
#'@@import combinat
@}

\section{Exchangeable multinomial model}\label{S:MLE} 

\subsection{Definitions} Let  $\mathbf{R}=(R_{1},\ldots, R_{K})^{T}$ follow an exchangeable multinomial distribution with $K+1$ categories.
We parameterize it by

\begin{equation} 
 \taurn{k}{n} = \Prob\big[\X_{\{1,\ldots,r_1\}}(O_1),  \ldots,
\X_{\{\sum_{i=1}^{k-1}r_i+1,\ldots,\sum_{i=1}^{k}r_i\}}(O_k)\big] \quad (k=1,\ldots,K),
\end{equation}
where $r_{i}\geq 0$ and $r_{1}+\cdots +r_k\leq n$. For notational convenience,
also let $\tau_{0,\ldots,0}=1$.

\subsection{Estimation} 
Consider $\taurn{K}{n}$ and its unconditional counterpart 
\begin{equation*} 
 \thetar{K}=\Prob\big[\X_{\{1,\ldots,r_1\}}(O_1), \ldots,
\X_{\{\sum_{i=1}^{K-1}r_i+1,\ldots,\sum_{i=1}^{K}r_i\}}(O_K)\big] = \sum_{n=\sum r_i}^{C}\taurn{K}{n} \Prob(N=n). 
\end{equation*}

If $\Arn{K}{n}$ denotes the number of clusters of size $n$ with response vector $(r_1,\ldots, r_K)$, then their non-parametric estimates are 
\begin{equation}  \label{E:mle} 
 \htaurn{K}{n}=\sum_{s_1,\ldots,s_{K}}\frac{\dbinom{n-\tsum{r_i}}{s_1,\ldots,s_{K}}}%
  {\dbinom{n}{r_1+s_1,\ldots,r_K+s_{K}}}\frac{\Arn{K}{n}}{M_{n}},% 
\end{equation}% 
and 
\begin{equation}  \label{E:thetahat} 
\hthetar{K}=\sum_{n=1}^M\sum_{s_1,\ldots,s_{K}}\frac{\dbinom{n-\tsum{r_i}}{s_1,\ldots,s_{K}}}%
  {\dbinom{n}{r_1+s_1,\ldots,r_K+s_{K}}}\frac{A_{r_1+s_1,\ldots,r_K+s_K|n}}{M}.
\end{equation}%


The function \texttt{tau} creates a ``look-up table'' for the MLEs. It returns either a list by treatment group
of either $K+1$ or $K$ dimensional arrays, depending on whether cluster-size specific estimates ($\tau$'s) or 
averaged estimates ($\theta$'s) are requested. For the cluster-size specific estimates the first dimension is
the cluster size. The calculation of $\theta$'s is done separately for each dose level, and thus each dose 
level uses a different sample-size distribution for averaging.


@o ..\R\ExchMultinomial.R @{
#'@@rdname CorrBin-internal
@< Define function for multinomial coefficient @>
tau <- function(cmdata, type=c("averaged","cluster")){
  type <- match.arg(type)
  
  @< Extract info from cmdata into variables @>
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
      @< Calculate averaged thetas @>
    } else {
      Mn <- xtabs(Freq ~ factor(ClusterSize, levels=1:M), data=cm1) 
      @< Calculate cluster-specific taus @>
    }
    
    # append treatment-specific result to result list
    res.trt <- list(res.trt)
    names(res.trt) <- trt
    res <- c(res, res.trt) 
  }
  res
}
@| tau @}


@d Extract info from cmdata into variables @{
  nc <- attr(cmdata, "ncat")
  nrespvars <- paste("NResp", 1:nc, sep=".")
  M <- max(cmdata$ClusterSize)
@}

First, we define the MLE averaged over cluster sizes. The \texttt{Calculate averaged thetas} macro
creates a $K$-dimensional array of $\thetar{K}(d)$ values.  
The implementation is based on combining the two summations
of the definition into one using $n=\sum_{i=1}^K r_i + \sum_{i=1}^K s_i + s_{K+1}$:

\begin{multline} 
\hthetar{K}=\sum_{n=1}^M\sum_{s_1,\ldots,s_{K}}\frac{\dbinom{n-\tsum{r_i}}{s_1,\ldots,s_{K}}}%
  {\dbinom{n}{r_1+s_1,\ldots,r_K+s_{K}}}\frac{A_{r_1+s_1,\ldots,r_K+s_K|n}}{M} \\
  = \sum_{s_1,\ldots,s_{K+1}}\frac{\dbinom{\tsum{s_i}}{s_1,\ldots,s_{K}}}%
  {\dbinom{\tsum r_i + \tsum{s_i}}{r_1+s_1,\ldots,r_K+s_{K}}}\frac{A_{r_1+s_1,\ldots,r_K+s_K|\tsum r_i + \tsum{s_i}}}{M}.
\end{multline}%

@D Calculate averaged thetas @{
    
res.trt <- array(NA, dim=rep(M+1, nc-1))
dimnames(res.trt) <- rep.int(list(0:M), nc-1) 
names(dimnames(res.trt)) <- paste("R", 1:(nc-1), sep="")
# indices for possible values of r
@<Simplex with sums @(idx @, M @, nc-1 @, idxsum @)@>
#indices for possible values of s 
# (one more column than for r - ensures summation over all n's)
@<Simplex with sums @(sidx @, M @, nc @, sidxsum @)@>
for (i in 1:nrow(idx)){
  r <- idx[i,]
  s.idx <- which(sidxsum <= M-sum(r))
  lower.idx <- sidx[s.idx, , drop=FALSE]
  upper.idx <- lower.idx + rep(c(r,0), each=nrow(lower.idx))
  res.trt[rbind(r)+1] <- 
    sum(mctab[lower.idx+1] / mctab[upper.idx+1] * atab[upper.idx+1]) / Mn
}
@}

Next, we define the MLEs specific for each cluster size. The macro \texttt{Calculate cluster-specific taus}
creates a $K+1$ dimensional array, with the cluster size as the first dimension.

@D Calculate cluster-specific taus @{
res.trt <- array(NA, dim=c(M, rep(M+1, nc-1))) #first dimension is 'n'
dimnames(res.trt) <- c(list(1:M), rep.int(list(0:M), nc-1)) 
names(dimnames(res.trt)) <- c("N",paste("R", 1:(nc-1), sep=""))
for (n in which(Mn > 0)){
  # indices for possible values of r
  @<Simplex with sums @(idx @, n @, nc-1 @, idxsum @)@>
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
@}


\section{Marginal compatibility}

Under marginal compatibility,
\begin{equation}
\pi_{\rvec|n} = \sum_{\tvec \in \V_M} h(\rvec, \tvec, n) \pi_{\tvec|M},
\end{equation}
where $h(\rvec, \tvec, n)  = \binom{\tvec}{\rvec}\binom{M-\sum t_i}{n-\sum r_i} \big/ \binom{M}{n} = 
\prod_{i=1}^K \binom{t_i}{r_i}\binom{M-\sum t_i}{n-\sum r_i} \big/ \binom{M}{n}$ and
$\V_n=\{(v_1,\ldots,v_K)\in \mathbb{N}^K \mid v_i \geq 0, \sum v_i \leq n\}$ is a $K$-dimensional simplex lattice with maximum
sum $n$.

\subsection{Estimation}

The following code implements the EM-algorithm for estimating the probabilities
of response assuming marginal compatibility. Let $(\rvec_i, n_i)$, $i=1,\ldots N$ denote
the observed data for a given dose level, where $i$ iterates
through the clusters, $n_i$ is the cluster size and 
$\rvec_i = (r_1,\ldots,r_K)$ is the observed number of responses of each type.

\begin{equation}\label{F:EMupdate0}
 \pi_{\tvec|M}^{(t+1)} = \frac{1}{N} \sum_{i=1}^{N} h(\rvec_{i},\tvec,n_{i})
             \frac{\pi^{(t)}_{\tvec|M}}{\pi^{(t)}_{\rvec_{i}|n_{i}}},
\end{equation}

The \texttt{mc.est.CMData} function implements the \texttt{mc.est} S3 method for \texttt{CMData} objects, 
returning a data frame with all  $\pi^{(g)}_{\rvec|n}, n=1,\ldots, M$ probabilities. The `hard' work is done
by the \texttt{mc.est.raw} function, which returns a list of matrices with  $\pi^{(g)}_{\rvec|M}$ values.
@O ..\R\ExchMultinomial.R
@{
#'@@rdname mc.est
#'@@method mc.est CMData
#'@@S3method mc.est CMData
#'@@export
#'@@param eps numeric; EM iterations proceed until the sum of squared changes fall below \code{eps}  


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
    fin <- fin[fin[last.resp] >= 0,]  #remove impossible clusters
    fin[c("Trt","ClusterSize", resp.vars1, last.resp, "Prob")]
}
@| mc.est.CMData @}

First we write a help-function that calculates all the probabilities
$\pi_{\rvec|n}$ given the set of $\theta_\rvec=\pi_{\rvec|M}$. While there are a variety
of ways doing this, we use a recursive formula:
\begin{equation}
\pi_{\rvec|n}  = \sum_{i=1}^K \frac{r_i+1}{n+1}\pi_{\rvec+\dvec_i|n+1} + \frac{n-\sum_ir_i+1}{n+1}\pi_{\rvec|n+1},
\end{equation}
where $\dvec_i$ is the $i$th coordinate basis vector (i.e.\ all its elements are 0, except the $i$th, which is 1).

The input for \texttt{Marginals} is a $K$-dimensional array of $\pi_{\rvec|M}$, and the output is a $(K+1)$-dimensional
array with the values of $\pi_{\rvec|n}$, $n=1,\ldots,M$ with cluster size as the first dimension

@O ..\R\ExchMultinomial.R
@{
#'@@rdname CorrBin-internal
Marginals <- function(theta){
  K <- length(dim(theta))
  M <- dim(theta)[1]-1
  
  res <- array(0, dim=c(M, rep(M+1, K)))
  dimnames(res) <- c(N=list(1:M), dimnames(theta))
  
  # indices for possible values of r
  @<Simplex with sums @(idx @, M @, K+1 @, clustersize @)@>
  idx <- idx[ , -1, drop=FALSE]  #remove (K+1)st category
  
  @< Initialize for cluster size M @>
  for (cs in seq.int(M-1,1)){
    @< Calculate values for cluster size cs... @>
  }
  
  res
}
@| Marginals @}

The initialization just copies over the values from \texttt{theta} to the appropriate dimension. Note that when indexing
the arrays, a ``+1'' is necessary since \texttt{idx} is 0-based.
@d Initialize for cluster size M @{
  curridx <- idx[clustersize==M, ,drop=FALSE]
  res[cbind(M, curridx+1)] <- theta[curridx+1]
@}
The iterative step initializes with the last term (with $\pi_{\rvec|n+1}$) and loops over the basis vectors.
@d Calculate values for cluster size cs given values for size cs+1 @{
  curridx <- idx[clustersize==cs, , drop=FALSE]
  res[cbind(cs, curridx+1)] <- (cs+1- rowSums(curridx))/(cs+1) * res[cbind(cs+1, curridx+1)]
  for (j in 1:K){
    lookidx <- curridx
    lookidx[ ,j] <- lookidx[ ,j] + 1   #add 1 to the j-th coordinate
    res[cbind(cs, curridx+1)] <- res[cbind(cs, curridx+1)] + 
                                 lookidx[,j]/(cs+1) * res[cbind(cs+1, lookidx+1)]
  }  
@}

The actual EM iterations are performed in \texttt{mc.est.raw}. 

@O ..\R\ExchMultinomial.R
@{
#'@@rdname CorrBin-internal
mc.est.raw <- function(object, ...) UseMethod("mc.est.raw")

#'@@method mc.est.raw CMData
#'@@S3method mc.est.raw CMData
mc.est.raw.CMData <- function(object, eps=1E-6, ...){
  cmdata <- object
  @< Extract info from cmdata into variables @>
  
  # indices for possible values of r with clustersize = M
  @<Simplex with sums @(idx @, M @, nc-1 @, idxsum @)@>

  res <- list()
  for (trt in levels(cmdata$Trt)){
    cm1 <- cmdata[cmdata$Trt==trt,]
    if (nrow(cm1) > 0){
      # observed freq lookup table
      atab <- array(0, dim=rep(M+1, nc))
      a.idx <- data.matrix(cm1[,nrespvars])
      atab[a.idx + 1] <- atab[a.idx + 1] + cm1$Freq
      Mn <- sum(cm1$Freq)
      
      @< MC estimates for given dose group @>
      
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
}@| mc.est.raw@}

Within each dose group, the algorithm iterates until the sum of squared changes of the parameters is smaller
than the selected threshold \texttt{eps}.
@D MC estimates for given dose group @{
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
    
    @< Calculate res.new - the value of res.trt for next iteration @>
	
    sqerror <- sum((res.new[idx+1] - res.trt[idx+1])^2)
	res.trt <- res.new 
  }
@}

The update of the $\pi_{\tvec|M}$ is performed based on \eqref{F:EMupdate0} rewritten to combine
clusters of the same type:
\begin{equation}\label{F:EMupdate1}
 \pi_{\tvec|M}^{(t+1)} = \frac{1}{N} \sum_{(\rvec,n)}\frac{A_{\rvec,n}}{\pi^{(t)}_{\rvec|n}} 
                                    h(\rvec,\tvec,n)\pi^{(t)}_{\tvec|M},
\end{equation}
looping through each cluster type ($\rvec, n$), and updating all $\pi_{\tvec|M}$ values compatible
with this type. The compatible $\tvec$ vectors have $t_i\geq r_i$, so they can be written in the form
$\tvec = \rvec + \svec$, where $s_i\geq 0$ and $\sum s_i \leq M-\sum r_i$.

@D Calculate res.new - the value of res.trt for next iteration @{
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
@}


\subsection{Manipulating estimates}
It is helpful to have functions that can convert the marginally compatible estimates from the $\pi$-based
form obtained in the estimates to the $\tau$'s and to extract the variance-covariance matrix and the 
correlation parameters.

The \texttt{tau.from.pi} function takes a $K$-dimensional array of $\pi_\rvec$ values, and returns a $K$-dimensional
array of $\tau_\rvec$ values using
\begin{equation} 
 \tau_\rvec = \sum_{\svec}\frac{\binom{n-\sum r_i}{\svec}}%
  {\binom{n}{\rvec+\svec}}\pi_{\rvec+\svec}.% 
\end{equation}

@O ..\R\ExchMultinomial.R
@{ 
#'@@rdname CorrBin-internal
tau.from.pi <- function(pimat){
  K <- length(dim(pimat))
  n <- dim(pimat)[1] - 1
  res <- array(NA, dim=rep(n+1, K)) 
  dimnames(res) <- rep.int(list(0:n), K) 
  names(dimnames(res)) <- paste("R", 1:K, sep="")

  # multinomial lookup table
  mctab <- mChooseTable(n, K+1, log=FALSE)
  
  # indices for possible values of r
  @<Simplex with sums @(idx @, n @, K @, idxsum @)@>
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
@| tau.from.pi @}

The \texttt{p.from.tau} function function takes a $K$-dimensional array of $\tau_\rvec$ values, and returns a vector 
of marginal probabilities of success $\tau_{\dvec_i}$.

@o ..\R\ExchMultinomial.R
@{
#'@@rdname CorrBin-internal
p.from.tau <- function(taumat){
  K <- length(dim(taumat))
  idx <- diag(nrow=K)
  taumat[rbind(idx+1)]
}    
@} 

The function \texttt{mc.corr} calculates the within- and between-outcome correlation coefficients for the
exchangeable model. It takes a $K$-dimensional array of $\pi_\rvec$ values, and returns a 2-dimensional
matrix of $\phi_{ij}$, $i,j=1, \ldots,K$ values using
\begin{equation} 
 \phi_{ij} = 
   \begin{cases}
   \big[\tau_{(2\dvec_i)} - \tau_{(\dvec_i)}^2\big]\big/\big[\tau_{(\dvec_i)} (1-\tau_{(\dvec_i)})\big] & i=j\\
  -\big[\tau_{(\dvec_i+\dvec_j)} - \tau_{(\dvec_i)}\tau_{(\dvec_j)}\big]\big/\big[\tau_{(\dvec_i)} \tau_{(\dvec_j)}\big] & i\ne j,\\
   \end{cases} 
\end{equation}
where $\dvec_i=(0,\ldots,0,\overbrace{1}^i,0,\ldots,0)$.
@O ..\R\ExchMultinomial.R
@{
#'@@rdname CorrBin-internal
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

@| corr.from.pi @}

\subsection{Testing marginal compatibility}
The \texttt{mc.test.chisq} function implements a generalization of the Cochran-Armitage trend test
for correlated multinomial data to test for marginal compatibility. Note that it only tests that the marginal probability of 
response $p_i$ does not depend on the cluster size for any category.

First, we define the test statistic for one group, and then add the resulting $\chi^2_{K}$-distributed test statistics
over the $G$ groups for an overall $G\,K$ degree of freedom test.

As above, let $(\rvec_i, n_i)$, $i=1,\ldots N$ denote
the observed data for a given dose level, where $i$ iterates
through the clusters, $n_i$ is the cluster size and 
$\rvec_i = (r_{i1},\ldots,r_{iK})$ is the observed number of responses of each type. Define the raw trend statistic for response $j$
as
\begin{equation}\label{E:rawmcstat}
X_j = \sum_{i=1}^N r_{ij} (c_{n_i} - \bar{c}), \quad j=1,\ldots,K,
\end{equation}
where $c_n$ are the scores for the Cochran-Armitage test usually chosen as $c_n=n-(M+1)/2$, and 
$\bar{c}_g=\big(\sum_{i=1}^{N}n_{i}c_{n_{i}}\big) \big/ \big(\sum_{i=1}^{N}n_{i}\big) = 
\sum_{n=1}^M M_nnc_n / \sum_{n=1}^M n M_n$ is the weighted
average of the scores ($M_n$ is the number of clusters of size $n$).

The covariance of two of these test statistics is
\begin{equation}
\sigma_{jk}=\text{Cov}(X_j, X_k) = 
\begin{cases}
\dsum_{i=1}^N (c_{n_i}-\bar{c})^2 n_i p_{j|n} (1-p_{j|n})[1+(n_i-1)\phi_{jj|n_i}], & j=k;\\ 
-\dsum_{i=1}^N (c_{n_i}-\bar{c})^2 n_i p_{j|n} p_{k|n}[1+(n_i-1)\phi_{jk|n_i}], & j\ne k,\\ 
\end{cases}
\end{equation}
where $p_{j|n}=\tau_{\dvec_j|n}$ is the probability of event type $O_j$ in clusters of size $n$. Under the
null hypothesis of marginal compatibility, the dependence of $p_{j|n}$ and $\phi_{jk|n}$ on $n$ can be removed: 
\begin{equation}\label{E:mcstatcov}
\sigma_{jk}=
\begin{cases}
 p_{j}(1-p_{j})\dsum_{i=1}^N (c_{n_i}-\bar{c})^2 n_i[1+(n_i-1)\phi_{jj}] = 
   p_{j}(1-p_{j})\dsum_{n=1}^M n M_n (c_n-\bar{c})^2 [1+(n-1)\phi_{jj}], & j=k;\\ 
-p_{j}p_{k}\dsum_{i=1}^N (c_{n_i}-\bar{c})^2 n_i [1+(n_i-1)\phi_{jk}] =
  -p_{j}p_{k}\dsum_{n=1}^M n M_n (c_n-\bar{c})^2 [1+(n-1)\phi_{jk}], & j\ne k,\\ 
\end{cases}
\end{equation}

The combined test statistic for the given dose group $g$ is
\begin{equation}
T^2_g = X_g' \Sigma_g^{-1} X_g \sim \chi^2_K \text{ under }H_0,
\end{equation}
where $X_g' =  (X_{g1}, \ldots, X_{gK})$, and $\Sigma_g=(\sigma_{gjk})_{K\times K}$ is its variance-covariance matrix defined by
\eqref{E:rawmcstat} and \eqref{E:mcstatcov}.
The unknown values of $p_j$ and $\phi_{jk}$ will be replaced by their estimates under marginal compatibility
$\hat{\tau}_{g\dvec_j}$ and $\hat{\phi}_{gjk}$.

The final test statistic is an independent combination of the statistics for each dose group:
\begin{equation}
T^2=\sum_{g=1}^G T_g^2 \sim \chi^2_{G\,K} \text{ under }H_0.
\end{equation}

@O ..\R\ExchMultinomial.R
@{
#'@@rdname mc.test.chisq
#'@@method mc.test.chisq CMData
#'@@S3method mc.test.chisq CMData
#'@@export

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
@| mc.test.chisq @}    


\section{Support functions}

The \texttt{Simplex with sums} macro creates a matrix (parameter 1) with rows containing the coordinates of an
integer lattice within a $d$-dimensional (parameter 3) simplex of size $n$ (parameter 2). That is all $d$-dimensional
vectors with non-negative elements with sum not exceeding $n$ are listed. The actual sums are saved in a vector (parameter 4).
Since this is a parametrized macro, it will expand to code, so no actual function calls will be made by the program.
This should reduce copying of the potentially large matrices.

@d Simplex with sums @{
   @1 <- hcube(rep(@2+1, @3))-1
   @4 <- rowSums(@1)
   @1 <- @1[@4 <= @2, ,drop=FALSE]  #remove impossible indices
   @4 <- @4[@4 <= @2]
@}

The \texttt{mChoose} function calculates the multinomial coefficient $\binom{n}{r_1,\ldots,r_K}$. The lower
part of the expression is passed as a vector. If its values add up to less than $n$, an additional value
is added. The function is not vectorized.

@d Define function for multinomial coefficient @{
    mChoose <- function(n, rvec, log=FALSE){
      rlast <- n - sum(rvec)
      rveclong <- c(rvec, rlast)
      if (any(rveclong < 0)) return(0)
      
      res <- lgamma(n + 1) - sum(lgamma(rveclong + 1))
      if (log) res else exp(res)
    }
@| mChoose @}

The \texttt{mChooseTable} function creates a lookup table of the multinomial coefficients 
with the number of categories $k$ and $n=\max \sum r_i$ given. The results is a $k$-dimensional array, with element
\texttt{[r1,\ldots,rK]} corresponding to $\binom{\sum (r_i-1)}{r_1-1,\ldots,r_k-1}$ (because the array is 1-indexed, while
$r_i$ can go from 0). The values in the array with coordinate sum exceeding $n$ are missing.
 
@o ..\R\ExchMultinomial.R @{
#'@@rdname CorrBin-internal
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
@}
\section{Files}

@f

\section{Macros}

@m

\section{Identifiers}

@u

\end{document}
