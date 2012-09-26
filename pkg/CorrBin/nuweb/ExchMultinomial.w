\documentclass[reqno]{amsart}
\usepackage[margin=1in]{geometry}
\usepackage[colorlinks=true,linkcolor=blue]{hyperref}
\renewcommand{\NWtarget}[2]{\hypertarget{#1}{#2}}
\renewcommand{\NWlink}[2]{\hyperlink{#1}{#2}}

\providecommand{\tsum}{\textstyle\sum}
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

First we need to set up the C file so that it can access the R internals.
@o ..\src\ReprodMultiCalcs.c
@{
#include <stdlib.h>
#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <R_ext/Applic.h>
@}

We will be using object of \texttt{CMData} class, which is defined in \texttt{CMData.w}.

We will also need to load support libraries.

@o ..\R\ExchMultinomial.R @{
  library(combinat)
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
@< Define function for multinomial coefficient @>
tau <- function(cmdata, type=c("averaged","cluster")){
  type <- match.arg(type)
  
  @< Extract info from cmdata into variables @>
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
@}

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

The actual EM iterations are performed in \texttt{mc.est}. 

@O ..\R\ExchMultinomial.R
@{
mc.est <- function(cmdata, eps=1E-6){
  @< Extract info from cmdata into variables @>
  
  # indices for possible values of r with clustersize = M
  @<Simplex with sums @(idx @, M @, nc-1 @, idxsum @)@>

  res <- list()
  for (trt in levels(cmdata$Trt)){
    cm1 <- subset(cmdata, Trt==trt)
    # observed freq lookup table
    atab <- array(0, dim=rep(M+1, nc))
    a.idx <- data.matrix(cm1[,nrespvars])
    atab[a.idx + 1] <- atab[a.idx + 1] + cm1$Freq
    Mn <- sum(cm1$Freq)
    
    @< MC estimates for given dose group @>
    
    # append treatment-specific result to result list
    res.trt <- list(res.trt)
    names(res.trt) <- trt
    res <- c(res, res.trt) 
  }
  res
}@| mc.est@}

Within each dose group, the algorithm iterates until the sum of squared changes of the parameters is smaller
than the selected threshold \texttt{eps}.
@D MC estimates for given dose group @{
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




\subsection{Testing marginal compatibility}
The \texttt{reprodmc.test.chisq} function implements Pang and Kuk's version of
the test for marginal compatibility. Note that it only tests that the marginal probability of 
response $p_i$ does not depend on the cluster size. The original test was only
defined for one group and the test statistic was compared to $\chi^2_1$ (or more
precisely, it was a z-test), however the test is easily generalized by adding 
the test statistics for the $G$ separate groups and using a $\chi^2_G$ distribution.

\begin{equation}
Z_g = \Big[\sum_{i=1}^{N_g} (c_{n_{g,i}} - \bar{c}_g) r_{g,i}\Big] \bigg/
  \Big[\hat{p}_g(1-\hat{p}_g)\sum_{i=1}^{N_g}n_{g,i}(c_{n_{g,i}} - \bar{c}_g)^2 \{1+(n_{g,i}-1)\hat{\rho}_g\}\Big]^{1/2},
\end{equation}
where $c_n$ are the scores for the Cochran-Armitage test usually chosen as $c_n=n-(M+1)/2$, 
$\bar{c}_g=\big(\sum_{i=1}^{N_g}n_{g,i}c_{n_{g,i}}\big) \big/ \big(\sum_{i=1}^{N_g}n_{g,i}\big)$ is a weighted
average of the scores; $\hat{p}_g=\big(\sum_{i=1}^{N_g}r_{g,i}\big) \big/ \big(\sum_{i=1}^{N_g}n_{g,i}\big)$ 
is the raw response probability, and 
$\hat{\rho}_g=1-\big[\sum_{i=1}^{N_g}(n_{g,i}-r_{g,i})r_{g,i}/n_{g,i}\big] \big/ 
\big[\hat{p}_g(1-\hat{p}_g)\sum_{i=1}^{N_g}(n_{g,i}-1)\big]$ is the Fleiss-Cuzack estimate of the intra-cluster
correlation for the $g$th treatment group. 

\begin{equation}
X^2=\sum_{g=1}^G Z_g^2 \sim \chi^2_G \text{ under }H_0.
\end{equation}

@O ../R/Reprod.R
@{
  
mc.test.chisq <- function(cbdata){
  cbdata <- subset(cbdata, Freq>0)
 
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
