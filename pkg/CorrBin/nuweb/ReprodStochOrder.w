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
\providecommand{\tsum}{\textstyle\sum}

\title{Testing for stochastic ordering assuming marginal compatibility}
\author{Aniko Szabo}
\date{\today}

\begin{document}
\maketitle


\section{Preliminaries}

First we need to set up the C file so that it can access the R internals.
@o ..\src\ReprodCalcs.c
@{
#include <stdlib.h>
#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <R_ext/Applic.h>
@}

We will be using object of \texttt{CBData} class, which is defined in \texttt{CBData.w}.


\section{Marginal compatibility}

\subsection{Estimation}

The following C code implements the EM-algorithm for estimating the probabilities
of response assuming marginal compatibility.

\begin{equation}
 \pi_{rM}^{(t+1)} = \frac{1}{N}
\sum_{i=1}^{N} h(r_{i},r,n_{i})
             \frac{\pi^{(t)}_{r,M}}{\pi^{(t)}_{r_{i},n_{i}}},
\end{equation}

First we write a help-function that calculates all the probabilities
$\pi_{r,n}$ given the set of $\theta_r=\pi_{r,M}$. While there are a variety
of ways doing this, we use a recursive formula:
\begin{equation}
\pi_{r,n}  = \frac{n+1-r}{n+1}\pi_{r,n+1} + \frac{r+1}{n+1}\pi_{r+1,n+1}
\end{equation}

@O ..\src\ReprodCalcs.c
@{
double static **Marginals(double* theta, int maxsize)
{
 double **res;
 int r, n;
 
 res = malloc((maxsize+1)*sizeof(double *));
 for (n=0; n<=maxsize; n++){
	 res[n] = calloc(n+1, sizeof(double));
 }
 
 for (r=0; r<=maxsize; r++){
	 res[maxsize][r] = theta[r];
 } 
 
 for (n=maxsize-1; n>=1; n--){
	 for (r=0; r<=n; r++){
		 res[n][r] = (n+1.0-r)/(n+1.0)*res[n+1][r] + (r+1.0)/(n+1.0)*res[n+1][r+1];
	 }
 }
 
 return res;	
}
@}

The actual EM iterations are performed in \texttt{ReprodEstimates}.

@O ..\src\ReprodCalcs.c
@{
SEXP ReprodEstimates(SEXP nvec, SEXP rvec, SEXP freqvec)
{
 double *theta, *thetanew, sqerror, **marg;
 int i, maxsize, nr, ntot, r, ri, ni, n, fri, start;
 SEXP res;
 const double eps=1e-8;
 
 nr = LENGTH(nvec);
 maxsize = 0;
 ntot = 0;
 for (i=0; i<nr; i++){
	 if (INTEGER(nvec)[i]>maxsize){
		 maxsize = INTEGER(nvec)[i];
	 }
	 ntot += INTEGER(freqvec)[i];
 }
		 
 theta = malloc((maxsize+1) * sizeof(double));
 thetanew = malloc((maxsize+1) * sizeof(double));
 //starting values
 for (r=0; r<=maxsize; r++){
	 theta[r] = 1.0/(maxsize+1);
 }
 sqerror = 1;
 //EM update
 while (sqerror>eps){
	 sqerror = 0;
	 marg = Marginals(theta, maxsize);
	 for (r=0; r<=maxsize; r++) thetanew[r] = 0;
	 for (i=0; i<nr; i++){
		 ri = INTEGER(rvec)[i];
		 ni = INTEGER(nvec)[i];
		 fri = INTEGER(freqvec)[i];
		 for (r=ri; r<=maxsize-ni+ri; r++){
			 thetanew[r] += choose(ni,ri)*choose(maxsize-ni,r-ri)*theta[r]*fri*1.0/
			                marg[ni][ri];
		 }
	 }
	 for (r=0; r<=maxsize; r++){
		 thetanew[r] = thetanew[r]/(ntot*choose(maxsize,r)*1.0);
		 sqerror += R_pow_di(thetanew[r]-theta[r],2);
		 theta[r] = thetanew[r];
	 }
	 for(n = 0; n <= maxsize; n++) free(marg[n]);
     free(marg);
 }

 PROTECT(res = allocMatrix(REALSXP, maxsize+1, maxsize));
 marg = Marginals(theta, maxsize);
 for (n=1, start=0; n<=maxsize; n++, start+=(maxsize+1)){
	 for (r=0; r<=n; r++){
		 REAL(res)[start+r] = marg[n][r];
	 }
	 for (r=n+1; r<=maxsize; r++){
		 REAL(res)[start+r] = NA_REAL;
	 }
 }
 
 for(n = 0; n <= maxsize; n++) free(marg[n]);
 free(marg);
 free(theta);
 free(thetanew);
 
 UNPROTECT(1);
 return res;
 	
}
@| ReprodEstimates @}

We will call the C function from R to calculate estimates separately by 
treatment group. For the package the compiled library needs to be loaded.

%In order to use the C code, it has to be first compiled into a
%DLL by calling \texttt{Rcmd SHLIB ReprodCalcs.c} from console. Then we load it
%into the R code\footnote{this is for a Windows system; on Linux/Unix a .so
%file would be created and that's what would be loaded}:

@o ../R/zzz.R
@{
.First.lib <- function(libname, pkgname){
  library.dynam("CorrBin", pkgname, libname)
}

.Last.lib <- function(){
  library.dynam.unload("CorrBin")
}
@}


@o ../R/Reprod.R
@{
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
@| mc.est@}

\subsection{Testing marginal compatibility}
The \texttt{reprodmc.test.chisq} function implements Pang and Kuk's version of
the test for marginal compatibility. Note that it only tests that the marginal probability of 
response $p_i$ does not depend on the cluster size. The original test was only
defined for one group and the test statistic was compared to $\chi^2_1$ (or more
precisely, it was a z-test), however the test is easily generalized by adding 
the test statistics for the $G$ separate groups and using a $\chi^2_G$ distribution.

@O ../R/Reprod.R
@{
  
mc.test.chisq <- function(cbdata){
  cbdata <- subset(cbdata, Freq>0)
 
  get.T <- function(x){
      max.size <- max(x$ClusterSize)
      K <- sum(x$Freq)  
      K.r <- xtabs(Freq~ClusterSize, data=x)
      scores <- -(max.size - (2*(1:max.size)-1))/2
      #n.t <- expected.n.unobserved(x)
      #p.hat <- sum((0:max.size)*n.t)/(max.size*K)
      a.r <- by(x, x$ClusterSize, function(z){sum(z$NResp*z$Freq)})[1:length(K.r)]
      cl.sizes <- as.numeric(names(a.r))
      N <- sum(cl.sizes*K.r)
      p.hat <- sum(a.r)/N
      p.hat.r <- a.r/(cl.sizes * K.r)
      T.stat <-  sum(cl.sizes * K.r * scores[cl.sizes]*p.hat.r)
      E.T.stat <- p.hat * sum(cl.sizes * K.r * scores[cl.sizes])
      rho.hat <- 1-sum(x$Freq*(x$ClusterSize-x$NResp)*x$NResp/x$ClusterSize)/((N-K)*p.hat*(1-p.hat))
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
   chi.list <- list(chi.sq=chis, p=pchisq(chis, df=1, lower.tail=FALSE))
   overall.chi <- sum(chis)
   overall.df <- length(chis)
   list(overall.chi=overall.chi, overall.p=pchisq(overall.chi, df=overall.df, lower.tail=FALSE), 
        individual=chi.list)
}
@| mc.test.chisq @}    


\section{Estimation under stochastic ordering}

We implement the estimation of the stochastically ordered MLE based on the
mixture representation following Hoff. The parameter \texttt{turn} controls the
extension to unimodal down-then-up orderings. The default value of \texttt{turn=1}
corresponds to the stochastically ordered $H_a$.

The actual calculations will be done in C;
the following function sets up the call from R. The input dataset should have
variables \texttt{NResp}, \texttt{ClusterSize},\texttt{Trt} and \texttt{Freq} -- 
it is most easily achieved by creating a ``CBData'' object.

The computational details are controlled by setting the \texttt{control} argument.
It should be a list with parameter settings; the simplest way to generate a correct
list is a call to the \texttt{mixControl()} function.
\texttt{S} gives the basis of the mixing distribution -- its rows are all the
possible non-decreasing vectors (see below for how it is obtained). \texttt{Q}
is the mixing distribution, it is a $G$-dimensional matrix ($G$ is the number 
of treatment groups). For the ISDM method it is initialized so that the constant vectors of $S$
are given equal probability ($H_0$ puts weight only on the constant vectors),
that is $Q[0,0,\ldots,0]= Q[1,1,\ldots,1]= \cdots =Q[N,N,\ldots,N]=1/(N+1)$. For
the EM method we have to put some weight on all possible $\bv$ vectors, because
it changes $Q$ only multiplicatively.
    
@O ../R/Reprod.R
@{

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
			           as.integer(verbose),PACKAGE="CorrBin"),
      ISDM = .Call("ReprodISDM", Q, S, tab, as.integer(max.iter), as.integer(max.directions),
                   as.double(eps),  as.integer(verbose),PACKAGE="CorrBin"))
 
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
@|mix.mc.mle @}

The values supplied in the call to \texttt{mixControl} replace the defaults and 
a list with all possible arguments is returned. The returned list is used as 
the control argument to the \texttt{mix.mc.mle} function.

The \texttt{method} argument allows to select either the EM, or the ISDM method.
@O ../R/Reprod.R @{
mixControl <- function(method=c("ISDM","EM"), eps=0.005, max.iter=5000, 
      max.directions=0, start=ifelse(method=="ISDM", "H0", "uniform"), verbose=FALSE){
  method <- match.arg(method)
  start <- match.arg(start, c("uniform","H0"))
  list(method = match.arg(method), eps = eps, max.iter = max.iter,
       max.directions = max.directions, start=start, verbose = verbose)
}
@| mixControl @}


The \texttt{makeSmatrix} function creates a matrix the rows of which are all the
possible non-decreasing vectors. The idea is the following: there is a one-to-one
correspondence between a non-decreasing sequence of length $k$ using the values
$0,1,\dots, n$ and combinations of $k$ elements out of $n+k$. Imagine $n$ balls and
$k$ sticks; each of the possible arrangements of theses $n+k$ objects into a row 
is a combination. We can transform it to a nondecreasing sequence by counting the
number of balls to the left of the first stick, second stick, etc. For example,
with $n=5$ and $k=4$, a possible arrangement is $\circ | \circ ||\circ\circ\circ|$; 
this corresponds to the sequence (1,2,2,5). The \texttt{GEN} function generates
all possible combinations of $k$ elements out of $N=n+k$ elements by recursion
on $N$, while \texttt{InsertComb} transforms the resulting combination to a non-decreasing
sequence and puts it into the output matrix. The code for generating the combinations
was written by  Joe Sawada, 1997  and obtained from the Combinatorial Object
Server (http://theory.cs.uvic.ca/inf/comb/CombinationsInfo.html).

@O ..\src\ReprodCalcs.c
@{
SEXP makeSmatrix(SEXP size, SEXP ntrt){
int *a, i, pos, nn, kk, nS;
SEXP res;
	
	void InsertComb (void){
		int i, val, step;
		val = 0;
		step = 0;
		for (i=1; i<=nn; i++) {
			if (a[i]==1){
				INTEGER(res)[pos+step]=val;
				step += nS;
			}
			else { //a[i]=0;
				val++;
			}
		}
	    pos++;
    }	    	
	
	void Comb(int j, int m) {
		if (j > nn) InsertComb();
		else {
			if (kk-m < nn-j+1) {
				a[j] = 0; Comb(j+1, m);
			}
			if (m<kk) {
				a[j] = 1; Comb(j+1, m+1);
			}
		}
	}

	
   nn = asInteger(size) + asInteger(ntrt);
   kk = asInteger(ntrt);
   a = calloc(nn+1, sizeof(int));
   for (i=1; i<=kk; i++) a[i]=1;

   nS = choose(nn,kk);
   PROTECT(res = allocMatrix(INTSXP, nS, kk));
   pos = 0;
   
   Comb(1, 0);
  
   UNPROTECT(1);
   free(a);
 
  return res;
}
@| makeSmatrix @}

\texttt{DownUpMatrix} extends the capabilities of \texttt{makeSmatrix} to generate
vectors that have $q$ non-increasing elements followed $k-q$ non-decreasing
elements.

The general outline of the algorithm is as follows ($n$ is \texttt{size}, 
$k$ is \texttt{ntrt}, and $q$ is \texttt{turn}):
@O ../R/Reprod.R
@{
DownUpMatrix <- function(size, ntrt, turn){
  if ((turn<1)|(turn>ntrt)) stop("turn should be between 1 and ntrt")
  @<Take care of turn=1 and turn=ntrt@>
  @<Generate non-increasing sequences of length turn with values $\leq$ size@>
  res2list <- list()
  for (sq in 0:size){
    @<Generate non-decreasing sequences of length ntrt-turn with values between sq and size@>
  }
  @<Combine the two parts of the sequences@>
}
@| DownUpMatrix @}

Non-increasing sequences are generated by subtracting non-decreasing sequences
from  $n$:
@d Generate non-increasing...@{
  res1 <- .Call("makeSmatrix", as.integer(size), as.integer(turn),PACKAGE="CorrBin")
  res1 <- size - res1;
@}

\texttt{res2list} will be a list of matrices defining non-decreasing sequences
for all possible starting points $s_q$. We get them by generating all non-
decreasing sequences with values $\leq n-s_q$ and then adding $s_q$. 
@d Generate non-decreasing sequences of length ntrt-turn ...@{
  S <- .Call("makeSmatrix", as.integer(size-sq), as.integer(ntrt-turn),PACKAGE="CorrBin")
  res2list <- c(res2list, list(sq+S))
@}

To match all possible starts with all possible ends, we split \texttt{res1} into
a list based on the last value ($s_q$) and do a Cartesian product (everything-to-
everything merge) with the corresponding elements of \texttt{res2list}:
@d Combine the two parts... @{
  res1list <- by(res1, res1[,turn], function(x)x)
  res <- mapply(merge, res1list, res2list, MoreArgs=list(by=NULL), SIMPLIFY=FALSE)
  res <- data.matrix(do.call(rbind, res))
  rownames(res) <- NULL
  colnames(res) <- NULL
  res
@}

When $q=1$ or $q=k$, no merging is needed and we can speed up the algorithm by
treating them as special cases:
@d Take care of turn=1...@{
  if (turn==1){
     res <- .Call("makeSmatrix", as.integer(size), as.integer(ntrt),PACKAGE="CorrBin")
     return(res)
  }
  if (turn==ntrt){
     res <- .Call("makeSmatrix", as.integer(size), as.integer(ntrt),PACKAGE="CorrBin")
     return(size - res)
  }
@}

\subsection{ISDM estimation}

In the ISDM algorithm, at each step we move toward a new mixing distribution
by maximizing
\begin{equation}
l\big(\alpha_0Q^{(t)}+\sum_{i}\alpha_i\delta_{\bv^i_t}\big),
\end{equation}
where $\bv^1_t, \ldots, \bv^m_t$ corresponding to the at
most $m_0$ largest positive values of the directional derivative 
$D_{Q^{(t)}}(\bv)$, where
\begin{equation}\label{E:D}
D_Q(\bv ) = \sum_{g,i}\frac{h(r_{g,i},v_g,n_{g,i})}{\sum_\bq
h(r_{g,i},q_g,n_{g,i})Q(\bq)}-N.
\end{equation}

First we define a variety of help-functions.

\texttt{GetTabElem} allows easy
access using three indices to the 3-dimensional data matrix \texttt{tab} that gets 
converted to a long vector when it is passed from R. Another option would have
been to set up a C-style array of pointers to pointers to rows.

@o ..\src\ReprodCalcs.c 
@{
double GetTabElem(SEXP tab, int size, int n, int r, int j){
	 return REAL(tab)[(n-1)+size*(r+(size+1)*j)];
 };
@}

\texttt{HyperTable} makes a 3-dimensional C-style array with the hyper-geometric probabilities
$h(i,j,k) =  {\binom{j}{i}}{\binom{N-j}{k-i}} / {\binom{N}{k}}$. R's \texttt{dhyper}
function does the actual calculations.

@O ..\src\ReprodCalcs.c
@{
double ***HyperTable(int size){
// dhyper(i, j, size-j, k), i=0:size; j=0:size; k=0:size	
 double ***res;
 int i, j, k;
 
 res = malloc((size+1)*sizeof(double*));
 for (i=0; i<=size; i++){
	 res[i] = malloc((size+1)*sizeof(double*));
	 for (j=0; j<=size; j++) res[i][j] = calloc(size+1, sizeof(double));
 }
 
 for (i=0; i<=size; i++){
	 for (j=i; j<=size; j++){
		 for (k=i; k<=size-j+i; k++)
			 res[i][j][k] = dhyper(i, j, size-j, k, 0); 
		 }
	 }

 return res;	
}
@| HyperTable @}


The matrix $S$ created by \texttt{makeSmatrix} contains the numeric values of the
non-decreasing vectors. \texttt{IndexVector} and \texttt{IndexVectorC} create a vector 
of indices corresponding to the location of each row of $S$ in a vectorized 
G-dimensional $N\times N\times\cdots\times N$ array ($Q$). That is if a row
of $S$ is $(0,1,4)$, then the corresponding element of the index vector will 
point to the [0,1,4]th element of a $N\times N\times N$ matrix that is represented
by a long (length$=N^3$) vector. This setup is needed as the number of dimensions
of $Q$ depends on the number of treatments, so it is not known at compile time.

The difference between the two functions is that
\texttt{IndexVector} works for an R-style matrix (that has been converted to a
vector during the transfer), while \texttt{IndexVectorC} works for a C-style matrix 
(an array of pointers to arrays).

@O ..\src\ReprodCalcs.c
@{
int static *IndexVectorC(int **S, int N, int G, int nrowS){
 int *idx, i, j;
 
 idx = calloc(nrowS, sizeof(int));
 
 for (j=G-1; j>=0; j--){
   for (i=0; i<nrowS; i++){
	   idx[i] = N*idx[i] + S[i][j];
   }
 }
 
 return idx;
}
int static *IndexVector(SEXP S, int N, int G, int nrowS){
 int *idx, i, j,  start;
 
 idx = calloc(nrowS, sizeof(int));
 
 for (j=G-1, start=(G-1)* nrowS; j>=0; j--, start -= nrowS){
   for (i=0; i<nrowS; i++){
	   idx[i] = N*idx[i] + INTEGER(S)[start+i];
   }
 }
 
 return idx;
}
@| IndexVector, IndexVectorC @}


\texttt{CalcMarginals} calculates the marginal probabilities in the denominator
of \eqref{E:Qupdate}
\begin{equation}
 Q_{gn}(r)= P(r \text{ responses } | \text{treatment }g, \text{cluster size }n)=
 \sum_\bq h(r,q_g,n) Q^{(t)}(\bq)
\end{equation}
for $j=1,\ldots, G$, $n=1,\ldots, N$, $r=0,\ldots,n$. It returns a 3-dimensional 
C-style matrix.

@O ..\src\ReprodCalcs.c 
@{
double ***CalcMarginals(SEXP S, SEXP Q, double ***ht, int *idx, int ntrt, int size, int nS){
   int j, i, n, x, start, sj;
	 double ***marg;
	 //ht is from HyperTable, it is passed to avoid recalculation
	 
	 marg = malloc(ntrt*sizeof(double*));
	 for (j=0; j<ntrt; j++){
		 marg[j] = malloc((size+1)*sizeof(double*));
		 for (n=1; n<=size; n++) marg[j][n] = calloc(n+1, sizeof(double));
	 }
	 
	 for (i=0; i<nS; i++){
		 for (j=0, start=0; j<ntrt; j++, start+=nS){
			 sj = INTEGER(S)[start+i];
			 for (n=1; n<=size; n++){
				 for (x=imax2(0,sj-size+n); x<=imin2(sj,n); x++){
					 marg[j][n][x] += REAL(Q)[idx[i]]*ht[x][n][sj];
				 }
			 }
		 }
	 }
  return marg;
 }
@| CalcMarginals @}
    
\texttt{CalcD} calculates the directional derivative of the log-likelihood into
the direction of each possible $\bv$ (that is for each row of $S$) as defined
in \eqref{E:D}. Note that the denominator contains the marginal probabilities 
defined above. The function does not have a return value, the first parameter 
$D$ is modified (or rather, the value it points to is modified).

@O ..\src\ReprodCalcs.c
@{
void CalcD(SEXP D, SEXP S, SEXP tab, int *idx, double ***ht, double ***marg, int ntrt, 
           int nS, int size, int ntot){
   int j, i, n, x, sj, start, t;
	 //ntot = sum(tab) -- passing it avoids having to recalculate it
	 for (i=0; i<nS; i++){
		 REAL(D)[idx[i]] = -ntot;
		 for (j=0, start=0; j<ntrt; j++, start+=nS){
			 sj = INTEGER(S)[start+i];
			 for (n=1; n<=size; n++){
				 for (x=imax2(0,sj-size+n); x<=imin2(sj,n); x++){
					 t = GetTabElem(tab,size,n,x,j);
					 if (t>0) REAL(D)[idx[i]] += t*ht[x][n][sj]/marg[j][n][x];
				 }
			 }
		 }
	 }
   	 
 }
@| calcD @}

\texttt{maxD} finds the largest value of the directional derivative $D$ --
it will be used for error rate estimation, since
\begin{equation}\label{E:relerror}
  l(\hat{Q}|\mathbf{X}) - l(Q|\mathbf{X}) \leq  \max_\bv D_Q(\bv).
\end{equation}

@o ..\src\ReprodCalcs.c
@{
double maxD(SEXP D, int *idx, int nS){
   int i;
   double currmax, val;
   
   currmax = 0;
   for (i=0; i<nS; i++){
	   val = REAL(D)[idx[i]];
	   if (val > currmax) currmax = val; 
   }
   return currmax;
}
@|MaxD @}


\texttt{CalcTopD} gets the vectors from $S$ corresponding to the at most 
`limit' positive D values (that is it calculates $\bv^1_t, \ldots, \bv^m_t$).
The return value is a C-style matrix with each row a non-decreasing vector.
The number $m$ of the selected vectors is returned through the pointer 
\texttt{nselect}.

@O ..\src\ReprodCalcs.c
@{
int **CalcTopD(SEXP D, SEXP S, int *idx, int limit, int *nselect, int ntrt, int nS){
	 int **res, nmx, i, j, pos, start;
	 double *posDvec,  dcut;
	 
	   // count the number of non-negative elements in D
     nmx = 0;
     for (i=0; i<nS; i++){
	     if (REAL(D)[idx[i]] >= 0){
		     nmx++;
	     }
     }
     if (nmx == 0){
	     res = 0;
	     *nselect = 0;
	     return res;
     }
     
     if (nmx > limit){ //find the limit-th largest D
		 posDvec = malloc(nmx*sizeof(double));
		 pos = 0;
	     for (i=0; i<nS; i++){
		     if (REAL(D)[idx[i]] >= 0){
			     posDvec[pos] = -REAL(D)[idx[i]];  //negation is needed, because rPsort uses increasing order
			     pos++;
		     }
	     }
	     rPsort(posDvec, nmx, limit);
	     dcut = -posDvec[limit];  //the cutoff for determining the limit-th largest values
	     free(posDvec);
     }
     else dcut = 0;
	 
     nmx = imin2(limit, nmx);

	 res = (int**) Calloc(nmx, int*);
	 pos = 0;
     for (i=0; i<nS; i++){
	     if (pos >= nmx) break;
	     if (REAL(D)[idx[i]] < dcut) continue;
	     res[pos] =(int*) Calloc(ntrt, int);  //copy the ith row of S
	     for (j=0, start=0; j<ntrt; j++, start+=nS) res[pos][j] = INTEGER(S)[start+i];
	     pos++;
     }

     *nselect = nmx;
     return res;	 
 } 
@| CalcTopD @}

\texttt{NegLogLik} calculates the negative log-likelihood at a potential 
new mixing distribution as a function of $(\alpha_0,\alpha_1,\ldots,\alpha_m)$, where
$\sum_{i=0}^m\alpha_i=1$ and $\alpha_i\geq 0$. To remove the sum-to-one constraint, 
we reparameterize to $\text{par}_i=\gamma_i=\alpha_{i}/\alpha_{0}\geq 0$, $i=1,\ldots,m$.  
With this reparametrization, $\alpha_0=1/(1+\sum\gamma_i)$ and $\alpha_i=\alpha_0\gamma_i$, $i=1,\ldots,m$,
and the constraints are $\gamma_i\geq 0$. The potential new mixing distribution
is 
\begin{equation*}
 Q^\gamma = \big(Q + \tsum_i\gamma_i\delta_{\bv^i}\big)/(1+\tsum_i\gamma_i)
\end{equation*}

The corresponding log-likelihood is
\begin{equation}\label{E:newloglik}
l(Q^\gamma) =
\sum_{g,n,r} A_{gnr}\log Q^\gamma_{gn}(r) \\=
\sum_{g,n,r} A_{gnr}\log\big(Q_{gn}(r)+\tsum_i\gamma_ih(r,v^i_g,n)\big) 
  - N\log(1+\tsum_i\gamma_i) ,
\end{equation}
where $A_{gnr}$ is the observed number of clusters of size $n$ with $r$ responses
in treatment group $g$.


\texttt{lmS} stores the $\bv^1_t, \ldots, \bv^m_t$ vectors calculated by \texttt{CalcTopD}.
This function will be defined within the scope of the \texttt{ReprodISDM} function
(because passing parameters through the \texttt{*ex} pointer is inconvenient),
so all the variables \texttt{tab, marg, ntot}, etc are available.

@D Define log-likelihood and its derivative
@{
double NegLogLik(int npar, double *par, void *ex){
	 //par[j] = (alpha_(j+1)/alpha_0), j=0,...,nmax-1
   int j, n, r, i, sj, x;
   double res, sum;
   
   res = 0;
   
   for (j=0; j<ntrt; j++){
	   for (n=1; n<=size; n++){
		   for (r=0; r<=n; r++){
			   x = GetTabElem(tab,size,n,r,j);
			   if (x>0){
				   sum = marg[j][n][r];
				   for (i=0; i<npar; i++){
					   sj = lmS[i][j];
					   sum += par[i]*ht[r][n][sj];
				   }
				   res += x*log(sum); 	 
			   }
		   }
	   }
   }
   sum = 0;
   for (i=0; i<npar; i++) sum += par[i];
   res -= ntot*log1p(sum);  //log1p(sum)=log(1+sum)
   
   if (!R_FINITE(res)){
     res = 1e60;
   }
   
   return (-res); }
@| NegLogLik @}

\texttt{NegLogLikDeriv} calculates the gradient vector of the negative log-likelihood
function defined above.
\begin{equation}
\frac{\partial l}{\partial \gamma_u} = \sum_{g,n,r}A_{gnr}\frac{h(r,v^u_g,n)}%
{Q_{gn}(r)+\tsum_i\gamma_ih(r,v^i_g,n)} - \frac{N}{1+\tsum_i\gamma_i}
\end{equation}

@D Define log-likelihood... @{
void NegLogLikDeriv(int npar, double *par, double *gr, void *ex){
   int j, n, r, i, sj, x;
   double alpha0, sum, ***denom;
      
   //prepare the shared denominators
   denom = malloc(ntrt*sizeof(double*));
   for (j=0; j<ntrt; j++){
	   denom[j] = malloc((size+1)*sizeof(double*));
	   for (n=1; n<=size; n++) denom[j][n] = calloc(n+1, sizeof(double));
   }
   for (j=0; j<ntrt; j++){
	   for (n=1; n<=size; n++){
		   for (r=0; r<=n; r++){
			   sum = marg[j][n][r];
			   for (i=0; i<npar; i++){
				   sj = lmS[i][j];
				   sum += par[i]*ht[r][n][sj];
			   }
			   denom[j][n][r] = sum; 	 
		   }
	   }
   }
   
   alpha0 = 1;
   for (i=0; i<npar; i++)  alpha0 += par[i];
   alpha0 = 1.0/alpha0;

   //calc the gradients 
   for (i=0; i<npar; i++){
	   sum = -ntot*alpha0;
	   for (j=0; j<ntrt; j++){
		   for (n=1; n<=size; n++){
			   for (r=0; r<=n; r++){
				   x = GetTabElem(tab,size,n,r,j);
				   if (x>0){
					   sj = lmS[i][j];
					   sum += x*ht[r][n][sj]/denom[j][n][r];
				   }
			   }
		   }
	   }
	   gr[i] = -sum;
   }
   
   for (j=0; j<ntrt; j++){
	   for (n=1; n<=size; n++) free(denom[j][n]);
	   free(denom[j]);
   }
   free(denom); 
}
@}

\texttt{UpdateQ} performs the updating step of the ISDM algorithm once the optimal
$\alpha_i$ values have been found (the input parameter \texttt{g} has the reparameterized 
\texttt{par} values at the optimum): the mixing
distribution is moved toward the values with the largest directional derivative.
\begin{equation}
Q^{(t+1)} =  \big(Q + \tsum_i\gamma_i\delta_{\bv^i}\big)/(1+\tsum_i\gamma_i)
\end{equation}
\texttt{lmS\_idx} is the index vector corresponding to the location in $Q$ of the 
selected directions.

@o ..\src\ReprodCalcs.c
@{
void UpdateQ(SEXP Q, double *g, int nS, int nmax, int *idx, int *lmS_idx){
   double alpha0;
   int i;
   
   alpha0 = 1;
   for (i=0; i<nmax; i++) alpha0 += g[i];
   alpha0 = 1/alpha0;
   
   for (i=0; i<nS; i++) REAL(Q)[idx[i]] *= alpha0;
   for (i=0; i<nmax; i++) REAL(Q)[lmS_idx[i]] += alpha0 * g[i];
}
@|UpdateQ @}

\texttt{UpdateMarginals} updates the marginal distribution calculated by
\texttt{CalcMarginals} after each iterative step -- this is faster than 
recalculating them again. \texttt{lmS} contains the selected directions.
\begin{equation}
Q^\gamma_{gn}(r) = \big(Q_{gn}(r) + \tsum_i\gamma_if(r,v^i_g,n)\big) \frac{1}{1+\tsum_i\gamma_i}
\end{equation}

@O ..\src\ReprodCalcs.c
@{
 void UpdateMarginals(double ***marg, double *g, double ***ht, int **lmS,
                      int ntrt, int size, int nmax){
   double alpha0;
   int i, j, n, r, sj;
    
    alpha0 = 1;
	for (i=0; i<nmax; i++){
		alpha0 += g[i];
	}
	alpha0 = 1/alpha0;
	
	for (j=0; j<ntrt; j++){
		for (n=1; n<=size; n++){
			for (r=0; r<=n; r++){
				for (i=0; i<nmax; i++){
					sj = lmS[i][j];
					marg[j][n][r] += g[i] * ht[r][n][sj];
				}
				marg[j][n][r] *= alpha0;
			}
		}
	}	 
}
@| UpdateMarginals @}

Finally, \texttt{ReprodISDM} sets up all the required matrices, and runs the
ISDM iterations while monitoring the error. The maximization of the log-likelihood
is done by the built-in R function \texttt{lbfgsb} which implements the quasi-Newton 
method of Byrd95 that allows boundary constraints for optimization. It is
actually a minimization routine, so the negative likelihood defined above is used.
The \texttt{MaxDirection} parameter controls the number of directions which
are considered at each step of the ISDM algorithm. Values of 0 (or less) mean
setting it to the number of non-empty cells in the data.

@O ..\src\ReprodCalcs.c
@{

SEXP ReprodISDM(SEXP Q, SEXP S, SEXP tab, SEXP MaxIter, SEXP MaxDirections, 
                SEXP eps, SEXP verbose){
 SEXP dims, D,  res, margSXP, tmp;
 int i, j, n, r, size, ntrt, *idx, nS, niter, **lmS, nmax, fncount, grcount, fail, 
     *boundtype,  limit, *lmS_idx, nenforced;
 double ***ht, ***marg, ntot, rel_error, *gamma, *lower, *upper, NLLmin;
 char msg[60];

 @<Define log-likelihood... @>
  
 PROTECT(dims = GET_DIM(tab));
   size = INTEGER(dims)[0];
   ntrt = INTEGER(dims)[2];
 UNPROTECT(1);
 
 PROTECT(dims = GET_DIM(S));
   nS = INTEGER(dims)[0];
 UNPROTECT(1);

 
 PROTECT(D = duplicate(Q));
 for (i=0; i<length(Q); i++) REAL(D)[i] = 0;
 
 idx = IndexVector(S, size+1, ntrt, nS);
 ht = HyperTable(size);

 ntot=0;
 for (i=0; i<length(tab); i++) ntot += REAL(tab)[i];

 limit = INTEGER(coerceVector(MaxDirections, INTSXP))[0];
 if (limit <= 0){  //set it to the number of non-empty cells
   limit=0;
   for (i=0; i<length(tab); i++){
 	   if (REAL(tab)[i]>0) limit++;
    }
 }

 marg = CalcMarginals(S, Q, ht, idx, ntrt, size, nS);
 CalcD(D, S, tab, idx, ht, marg, ntrt, nS, size, ntot);
 
 rel_error = maxD(D, idx, nS);
 niter = 0;
 nenforced = 0;
 
 while (niter < asInteger(MaxIter) && rel_error > asReal(eps)){
	 R_CheckUserInterrupt();
	 niter++;
     
	 lmS = CalcTopD(D, S, idx, limit, &nmax, ntrt, nS);
	 lmS_idx = IndexVectorC(lmS, size+1, ntrt, nmax);
	 
	 if (nmax == limit) nenforced++;
	 
	 gamma = (double*) Calloc(nmax, double); 
	 lower = (double*) Calloc(nmax, double); 
	 upper = (double*) Calloc(nmax, double); 
	 boundtype = (int*) Calloc(nmax, int); 
	 
 	 for (i=0; i<nmax; i++){
		 gamma[i] = 0;
		 lower[i] = 0;
		 upper[i] = imin2(1e6/nmax, 100); // => alpha0>1e-6
		 boundtype[i] = 1; //lower  bound only
	 }
	 
	 
 	 lbfgsb(nmax, 5, gamma, lower, upper, boundtype, &NLLmin, NegLogLik, NegLogLikDeriv,
 	        &fail, 0, 1e5, 0, &fncount, &grcount, 1000, msg, asInteger(verbose), 10);
 	        
	 UpdateMarginals(marg, gamma, ht, lmS, ntrt, size, nmax);
   CalcD(D, S, tab, idx, ht, marg, ntrt, nS, size, ntot);
	 UpdateQ(Q, gamma, nS, nmax, idx, lmS_idx); //only needed to be able to return Q
	 rel_error = maxD(D, idx, nS);
	 if (asInteger(verbose)==1)
	   Rprintf("Step %d, rel.error=%f, NLL=%f\n", niter, rel_error, NLLmin);

	  	 
	 Free(gamma);
	 Free(lower);
	 Free(upper);
	 Free(boundtype);
	 for (i=0; i<nmax; i++){
		 Free(lmS[i]);
	 }
	 Free(lmS);
	 free(lmS_idx);
 }
 
 free(idx);
 for(j=0; j<=size; j++){
	  for (n=0; n<=size; n++) free(ht[j][n]);
	  free(ht[j]);
  }
 free(ht);
 
 PROTECT(margSXP = allocVector(REALSXP, ntrt*size*(size+1)));
 PROTECT(dims = allocVector(INTSXP,3));
    INTEGER(dims)[0] = size+1;     INTEGER(dims)[1] = size;     INTEGER(dims)[2] = ntrt;   
    i = 0;
    for(j=0; j<ntrt; j++){
	  for (n=1; n<=size; n++){ 
		  for (r=0; r<=n; r++){
			  REAL(margSXP)[i] = marg[j][n][r];
			  i++;
		  }
		  for (r=n+1; r<=size; r++){
			  REAL(margSXP)[i] = NA_REAL;
			  i++;
		  }
	  }
   }
 dimgets(margSXP, dims);
 UNPROTECT(1);    //dims
 for(j=0; j<ntrt; j++){
	  for (n=1; n<=size; n++) free(marg[j][n]);
	  free(marg[j]);
  }
 free(marg);

 
 PROTECT(res = allocVector(VECSXP,5));
  SET_VECTOR_ELT(res, 0, margSXP);
  SET_VECTOR_ELT(res, 1, Q);
  SET_VECTOR_ELT(res, 2, D);
  SET_VECTOR_ELT(res, 3, ScalarReal(-NLLmin));
  PROTECT(tmp = allocVector(REALSXP,2));
   REAL(tmp)[0] = rel_error;
   REAL(tmp)[1] = niter;
  SET_VECTOR_ELT(res, 4, tmp);
 UNPROTECT(4);   //tmp, res, margSXP, D

 if (asInteger(verbose)==1) 
   Rprintf("Limit=%d; %d Iterations; Limit enforced %d times (%4.2f percent)\n", 
           limit, niter, nenforced, nenforced*100.0/niter);

 return res;

} 

@| ReprodISDM @}

\subsection{EM estimation}
The \texttt{MixReprodQ} function implements the EM-based fitting:
\begin{equation}\label{E:Qupdate}
Q^{(t+1)}(\bv)=\frac{1}{N}\sum_{g,i}\frac{h(r_{g,i},v_g,n_{g,i}) Q^{(t)}(\bv)}{\sum_\bq
h(r_{g,i},q_g,n_{g,i}) Q^{(t)}(\bq)}.
\end{equation}
Note that from \eqref{E:D}
\begin{equation}
 Q^{(t+1)}(\bv) = Q^{(t)}(\bv) + \frac{1}{N}D_Q(\bv)Q^{(t)}(\bv),
\end{equation}
so we will be able to use the \texttt{CalcMarginals} and \texttt{CalcD} functions
defined earlier.

\texttt{UpdateReprodQ} performs the updating step defined in \eqref{E:Qupdate}.

@O ..\src\ReprodCalcs.c
@{
SEXP UpdateReprodQ(SEXP Q, SEXP S, SEXP tab, int size, int ntrt, int nS, 
                   double*** ht, int* idx){
 int ntot, i, j, n;
 double ***marg;
 SEXP resobj, D;
 
 PROTECT(resobj = duplicate(Q));
 for (i=0; i<length(Q); i++) REAL(resobj)[i] = 0;
 
 ntot=0;
 for (i=0; i<length(tab); i++) ntot += REAL(tab)[i];
 
 marg = CalcMarginals(S, Q, ht, idx, ntrt, size, nS);

 PROTECT(D = duplicate(Q));
 for (i=0; i<length(Q); i++) REAL(D)[i] = 0;
 CalcD(D, S, tab, idx, ht, marg, ntrt, nS, size, ntot);
  
// update Q values  
 for (i=0; i<length(Q); i++){
   REAL(resobj)[i] = REAL(Q)[i] * (1+REAL(D)[i]/ntot);
 }

 //cleanup
 for(j=0; j<ntrt; j++){
	  for (n=1; n<=size; n++) free(marg[j][n]);
	  free(marg[j]);
 }
 free(marg);

 
 UNPROTECT(2);
 return resobj;
}
@| UpdateReprodQ @}

And, finally, \texttt{MixReprodQ} performs the actual EM iterations. \texttt{eps}
controls the precision of the estimate of the log-likelihood. Since
\begin{equation}\label{E:relerror_EM}
  l(\hat{Q}|\mathbf{X}) - l(Q^{(t)}|\mathbf{X}) \leq N 
	     \frac{Q^{(t+1)}(\bv )-Q^{(t)}(\bv )}{Q^{(t)}(\bv )},
\end{equation}
we control the largest relative change in $Q$ during updating. If the \texttt{verbose}
option is selected, the current value of the relative change is displayed at
every 10th iteration.

@O ..\src\ReprodCalcs.c
@{
SEXP MixReprodQ(SEXP Q, SEXP S, SEXP tab, SEXP MaxIter, SEXP eps, SEXP verbose){
 double rel_error, ***ht, ntot, re, ***marg, loglik;
 int niter, size, ntrt, nS, i, j, n, r, *idx, Qlength;
 SEXP D, Qnew, resQ, dims, resobj, tmp, margSXP, Qdims;

 PROTECT(Qdims = GET_DIM(Q));
 Qlength = 1;
 for (i=0; i<length(Qdims); i++) Qlength *= INTEGER(Qdims)[i];
 
 PROTECT(resQ = duplicate(Q));
 
 PROTECT(dims = GET_DIM(tab));
   size = INTEGER(dims)[0];
   ntrt = INTEGER(dims)[2];
 UNPROTECT(1);  //dims
 
 PROTECT(dims = GET_DIM(S));
   nS = INTEGER(dims)[0];
 UNPROTECT(1); //dims
 
 ntot=0;
 for (i=0; i<length(tab); i++) ntot += REAL(tab)[i];

 idx = IndexVector(S, size+1, ntrt, nS);
 ht = HyperTable(size);

 PROTECT(tmp = allocVector(REALSXP, 2));
 rel_error = 1;
 niter = 0;
 while ((niter<asInteger(MaxIter))&&(rel_error>asReal(eps))){
	  R_CheckUserInterrupt();
		niter++;
    PROTECT(Qnew = UpdateReprodQ(resQ, S, tab, size, ntrt, nS, ht, idx));
    rel_error = 0;
    for (i=0; i<length(Qnew); i++){
	    if (REAL(resQ)[i]>0){
		    re = ntot*(REAL(Qnew)[i]-REAL(resQ)[i])/REAL(resQ)[i];
		    if (rel_error < re) rel_error = re;
	    }
	    REAL(resQ)[i] = REAL(Qnew)[i];
    }
    UNPROTECT(1); //Qnew
    if ((asInteger(verbose) == 1)&&(niter%10 == 1)){
	    REAL(tmp)[1] = rel_error;
	    REAL(tmp)[0] = niter;
        PrintValue(tmp);
     } 

  }	
 UNPROTECT(1);  //tmp	
 
 //calculate ML estimates for the output
 marg = CalcMarginals(S, resQ, ht, idx, ntrt, size, nS);
 PROTECT(margSXP = allocVector(REALSXP, ntrt*size*(size+1)));
 PROTECT(dims = allocVector(INTSXP,3));
    INTEGER(dims)[0] = size+1;     INTEGER(dims)[1] = size;     INTEGER(dims)[2] = ntrt;   
    i = 0;
    for(j=0; j<ntrt; j++){
	  for (n=1; n<=size; n++){ 
		  for (r=0; r<=n; r++){
			  REAL(margSXP)[i] = marg[j][n][r];
			  i++;
		  }
		  for (r=n+1; r<=size; r++){
			  REAL(margSXP)[i] = NA_REAL;
			  i++;
		  }
	  }
   }            
 dimgets(margSXP, dims);   
 UNPROTECT(1); //dims

 PROTECT(D = allocVector(REALSXP, Qlength));
 dimgets(D, Qdims);
 for (i=0; i<Qlength; i++) REAL(D)[i] = 0;
 CalcD(D, S, tab, idx, ht, marg, ntrt, nS, size, ntot); 

 //calculate log-likelihood for the output
 loglik = 0;
 for(j=0; j<ntrt; j++){
  for (n=1; n<=size; n++){ 
	  for (r=0; r<=n; r++){
	    loglik += GetTabElem(tab, size, n, r, j)*log(marg[j][n][r]);
     }
   }
  }

 for(j=0; j<ntrt; j++){
	  for (n=1; n<=size; n++) free(marg[j][n]);
	  free(marg[j]);
  }
 free(marg);


 for(j=0; j<=size; j++){
	  for (n=0; n<=size; n++) free(ht[j][n]);
	  free(ht[j]);
  }
 free(ht);
 free(idx);
 
 PROTECT(resobj = allocVector(VECSXP, 5));
   SET_VECTOR_ELT(resobj, 0, margSXP);
   SET_VECTOR_ELT(resobj, 1, resQ);
   SET_VECTOR_ELT(resobj, 2, D);
   SET_VECTOR_ELT(resobj, 3, ScalarReal(loglik));
	 PROTECT(tmp = allocVector(REALSXP,2));
	   REAL(tmp)[0] = rel_error;
	   REAL(tmp)[1] = niter;
	 SET_VECTOR_ELT(resobj, 4, tmp);
 UNPROTECT(6); //tmp, resobj,D, margSXP, resQ, Qdims
 
 return resobj;

}

@}

\section{Testing for stochastic ordering}

Now that we have the MLEs, we implement the likelihood-ratio test. 

\subsection{Testing against a global null}
First  we consider testing
\begin{equation*}
H_0: \bpi_1=\cdots =\bpi_G \quad\text{versus}\quad
H_a:\bpi_1 \leqst \cdots \leqst \bpi_G
\end{equation*}
We incorporate the \texttt{turn} parameter to allow fitting and testing umbrella orderings.

The \texttt{SO.LRT} computes the likelihood-ratio test statistic.
@O ../R/Reprod.R
@{
SO.LRT <- function(cbdata, control=mixControl()){
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
  res <- mix.mc.mle(cbdata, control=control)
  ll1 <- res$loglik
  lrt <- 2*(ll1 - ll0)
  attr(lrt, "ll0") <- ll0
  attr(lrt, "ll1") <- ll1
  lrt
 }
  
@}

To get the p-value for the LRT, we use a permutation testing approach. The \texttt{SO.trend.test} function
relies on the \texttt{boot} library for this, thus both the permutation test
is straightforward to implement. \texttt{R}
specifies the number of resamples, \texttt{method} could be either ``ISDM'' or 
``EM'', and \texttt{eps} is the precision of the LRT estimate.

@O ../R/Reprod.R
@{
SO.trend.test <- function(cbdata, R=100, control=mixControl()){
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
@| SO.trend.test @}




\subsection{Universal trend test function}
\texttt{trend.test} provides a common interface for the Rao-Scott, GEE, and stochastic
order based trend tests.

@O ../R/Reprod.R @{
trend.test <- function(cbdata, test=c("RS","GEE","GEEtrend","GEEall","SO"), exact=test=="SO", 
                       R=100, control=mixControl()){ 
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
@}

\subsection{Finding the NOAEL}

The NOSTASOT dose is the No-Statistical-Significance-Of-Trend dose -- the largest dose 
at which no trend in the rate of adverse events has been observed. 

@O ../R/Reprod.R @{
NOSTASOT <- function(cbdata, test=c("RS","GEE","GEEtrend","GEEall","SO"), exact=test=="SO",
                     R=100, sig.level=0.05, control=mixControl()){
   ntrt <- nlevels(cbdata$Trt)
   control.gr <- levels(cbdata$Trt)[1]
   p.vec <- array(NA, ntrt-1)
   names(p.vec) <- levels(cbdata$Trt)[-1]
   NOSTASOT.found <- FALSE
   curr.gr.idx <- ntrt
   curr.gr <- levels(cbdata$Trt)[ntrt]
   
   while (!NOSTASOT.found & (curr.gr.idx>1)){
     d1 <- subset(cbdata, unclass(Trt)<=curr.gr.idx)
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
@| NOSTASOT @}

\section{Files}

@f

\section{Macros}

@m

\section{Identifiers}

@u

\end{document}
