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

\title{Correlated multinomial data}
\author{Aniko Szabo}
\date{\today}


\begin{document}
\begin{abstract} We define a class for describing data from toxicology experiments with multinomial outcomes
and implement fitting of a variety of existing models and trend tests.
\end{abstract}
\maketitle

\section{Defining \texttt{CMData} -- a class for \textbf{C}lustered \textbf{M}ultinomial \textbf{Data}}
We start with defining an S3 class describing data from toxicology experiments with multinomial outcomes of type
$1, 2, \ldots, K+1$. Here $K$ denotes the ``degrees of freedom'' of the outcome. $K=1$ corresponds to binary data. The
class is a data frame with the following columns:

\begin{description}
\item[Trt] a factor defining (treatment) groups
\item[ClusterSize] an integer-valued variable defining the cluster size
\item[NResp.1--NResp.K] $K$ integer-valued variables defining the number of responses of type $1,2,\ldots,K$
\item[Freq]  an integer-valued  variable defining frequency for each
Trt/ClusterSize/NResp.1/$\cdots$/NResp.K combination
\end{description}

\texttt{CMData} converts a data frame to a CMData object. \texttt{x}
is the input data frame; \texttt{trt}, \texttt{clustersize},  and
\texttt{freq} could be strings or column indices defining the appropriate
variable in \texttt{x}. \texttt{nresp} should be a vector of variable names or column indices of length $K$.
\texttt{freq} can also be NULL, in which case it is 
assumed that each combination has frequency 1.
@O ../R/CMData.R @{
CMData <- function(x, trt, clustersize, nresp, freq=NULL){
  if (!is.data.frame(x)) stop("x has to be a data frame")
  nms <- names(x)
  K <- length(nresp)
  
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
  nrespvar <- sapply(nresp, process.var)
  if (K > 1) colnames(nrespvar) <- 1:K
  if (is.null(freq)) freqvar <- rep(1, nrow(x))
  else freqvar <- process.var(freq)
  
  d <- data.frame(Trt=trtvar, ClusterSize=csvar, NResp=nrespvar, Freq=freqvar)
  nrespnames <- grep("NResp", names(d), value=TRUE)
  d <- aggregate(d$Freq, d[,c("Trt", "ClusterSize", nrespnames)], sum)
  names(d)[length(names(d))] <- "Freq"
  
  attr(d, "K") <- K
  class(d) <- c("CMData", "data.frame")
  d}
@| CMdata.data.frame  @}


The \texttt{read.CMData} function reads in clustered multinomial data from a tab-delimited
text file. There are two basic data format options:  etiher the counts of responses of all categories are given (and the
cluster size is the sum of these counts), or  the total cluster size is given with the counts of all but one category.
The first column should always give the treatment group, then either the counts for each category (first option, chosen by setting 
\texttt{with.clustersize = FALSE}), or the size of the cluster followed by the counts for all but one category (second option,
chosen by setting \texttt{with.clustersize = TRUE}). Optionally, a last column could
give the number of times the given combination occurs in the data.

@o ../R/CMData.R
@{
read.CMData <- function(file, with.clustersize=TRUE, with.freq=TRUE, ...){
  d <- read.table(file, ...)
  K <- ncol(d) - with.freq - 2  #subtracting Trt & either ClusterSize or last category column
  nrespvars <- paste("NResp", 1:K, sep=".")
  if (with.clustersize) 
    names(d) <- c("Trt","ClusterSize", nrespvars, if (with.freq) "Freq"))
  else {
   names(d) <- c("Trt", nrespvars, "LastCat", if (with.freq) "Freq"))
   d$ClusterSize <- rowSums(d[, c(nrespvars, "LastCat")]
   d$LastCat <- NULL
  }
  d <- CMData(d, trt="Trt", clustersize="ClusterSize", nresp=nrespvars, 
              freq=if (with.freq) "Freq" else NULL)
  d}
@| read.CMData @}

\texttt{unwrap.CMData} is a utility function that reformats a CMData object so that
each row is one observation (instead of one cluster). A new `ID' variable is added
to indicate clusters. 

@O ../R/CMData.R @{

unwrap.CMData <- function(cmdata){
  #unwrap Freq variable
  freqs <- rep(1:nrow(cmdata), cmdata$Freq)
  cm1 <- cmdata[freqs,]
  cm1$Freq <- NULL
  
  #create ID
  cm1$ID <- factor(1:nrow(cm1))
  
  #create last category count
  nrespvars <- grep("NResp", names(cm1), value=TRUE)
  cm1$LastCat <- cm1$ClusterSize - rowSums(cm1[,nrespvars])
  K <- length(nrespvars)
  
  #reshape to have one row per category within cluster
  cm2 <- reshape(cm1, direction="long", varying=list(c(nrespvars, "LastCat")), v.names="Count",
                 idvar="ID", timevar="Outcome", times=1:(K+1))
  
  #unwrap categories
  counts <- rep(1:nrow(cm2), cm2$Count)
  res <- cm2[counts,]
  cm2$Count <- NULL

  res[order(res$ID),]
  }
@| unwrap.CMData @}

\end{document}
