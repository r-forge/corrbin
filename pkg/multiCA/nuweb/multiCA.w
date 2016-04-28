\documentclass[reqno]{amsart}
\usepackage[margin=1in]{geometry}
\usepackage[colorlinks=true,linkcolor=blue]{hyperref}
\renewcommand{\NWtarget}[2]{\hypertarget{#1}{#2}}
\renewcommand{\NWlink}[2]{\hyperlink{#1}{#2}}

\newcommand{\pvec}{\mathbf{p}}
\def\T{{ \mathrm{\scriptscriptstyle T} }}
\newcommand{\setJ}{\mathcal{J}}
\newtheorem{theorem}{Theorem}


\title{Test for trend with a multinomial outcome}
\author{Aniko Szabo}
\date{\today}

\begin{document}
\maketitle


\section{Introduction}
-
Consider a study in which a multinomial outcome with $K$ possible unordered values is measured in subjects belonging to one of $G$ ordered groups.
The size of each group, $n_{i\cdot}$, is defined by the study design, and will be treated as fixed.
Let $\pvec_i=(p_{i1},\ldots,p_{iK})^\T$ denote the probabilities of the multinomial outcomes in the $i$th group. The hypothesis of interest is to
evaluate the homogeneity of these probabilities across the groups with a targeted alternative of a trend in at least one of the categories.
Formally, we consider testing $H_0 = \bigcap_{j=1}^K H_{0j}$ versus $H_1 = \bigcup_{j=1}^K H_{1j}$, where
\begin{equation}
\begin{aligned}
  H_{0j}& :  p_{1j}=\cdots = p_{Gj} \\
  H_{1j}& : p_{1j} \leq \cdots \leq p_{Gj} \text{ or }  p_{1j} \geq \cdots \geq p_{Gj} \text{ with at least one inequality}
\end{aligned}
\end{equation}

The test is based on the following result:
\begin{theorem}\label{Th:partial}
Let $\setJ \subset \{1,\ldots,K\}$, then under $H_{0\setJ}=\bigcap_{j\in\setJ}H_{0j}$ as $N\rightarrow\infty$
\begin{equation}
  W_\setJ =  \sum_{j\in\setJ} (1-p_{\cdot j})T^2_j + \big(\sum_{j\in\setJ} p_{\cdot j}\big) T^2_\setJ \xrightarrow{d} \chi^2_{d},
\end{equation}
where $d = \min(|\setJ|, K-1)$,  $T_\setJ = [\sum_{i=1}^G \sum_{j\in\setJ}n_{ij}(c_i-\bar{c})] \big/ \sqrt{p_{\cdot \setJ}(1-p_{\cdot \setJ})s^2}$ denotes the Cochran-Armitage trend test statistic for testing for marginal trend in $p_{i\setJ}=\sum_{j\in\setJ}p_{ij}$, $i=1,\ldots, G$.
\end{theorem}

\section{Implementing the overall test}

The main \texttt{multiCA.test} function is a generic, with methods for a matrix and formula input.

@o ../R/aaa-generics.R
@{
#'Multinomial Cochran-Armitage trend test
#'
#'The \code{multiCA.test} performs a multinomial generalization of the 
#' Cochran-Armitage trend test.
#'
#'
#'@@export
#'@@param x a two-dimensional matrix or a formula
#'@@param \dots other arguments 
#'@@return an object of class "htest" with the results of the test
#'@@author Aniko Szabo
#'@@references Szabo, A. (2016) Test for trend with a multinomial outcome.  
#'@@keywords nonparametric 
#'@@examples
#'
#'data(stroke)
#'## using formula interface
#'multiCA.test(Type ~ Year, weights=Freq, data=stroke)
#'
#'## using matrix interface
#'strk.mat <- xtabs(Freq ~ Type + Year, data=stroke)
#'multiCA.test(strk.mat)
#'
#'@@name multiCA.test

multiCA.test <- function(x,...) UseMethod("multiCA.test")
 
@}

The default method uses a two-dimensional contingency matrix with the outcomes as rows and ordered groups as columns.
@O ../R/multiCA.R @{
#'@@rdname multiCA.test
#'@@method multiCA.test default
#'@@param scores numeric vector of the same length as the number of ordered groups. Defaults to linearly increasing values
#'@@param outcomes integer or character vector defining the set of outcomes (by row index or row name) over which the trend should be tested. Defaults to all outcomes.
#'@@export

multiCA.test.default <- function(x, scores=1:ncol(x), outcomes=1:nrow(x),
...){
  if (!is.matrix(x)) {
    cat(str(x))
    stop("x should be a two-dimensional matrix")
}
  if (length(scores) != ncol(x)) stop("The length of the score vector should equal the number of columns of x")

  K <- nrow(x)
  full <- length(outcomes) == K  #full test
  
  nidot <- apply(x, 2, sum)
  n <- sum(nidot)
  
  cbar <- sum(nidot * scores)/n
  
  s2 <- sum(nidot * (scores - cbar)^2)
  pdot <- prop.table(rowSums(x))[outcomes]
  nonz <- (pdot > 0)
  
  if (!any(nonz)) return(1)
  
  X <- x[outcomes, ,drop=FALSE] %*% (scores - cbar)
  
  if (full || sum(pdot) >= 1){
    Tt <- ( sum(X[nonz]^2 / pdot[nonz])) / s2
  } else {
    Tt <- (sum(X)^2 / (1-sum(pdot)) + sum(X[nonz]^2 / pdot[nonz])) / s2
  }
  names(Tt) <- "W"

  df <- length(outcomes) - full
  names(df) <- "df"
  p.value <- pchisq(Tt, df=df, lower.tail=FALSE)

  res <- list(statistic = Tt, parameter = df, p.value = p.value, 
              method="Multinomial Cochran-Armitage trend test",
              data.name = deparse(substitute(x)))
  class(res) <- "htest"
  return(res)  
}
@| multiCA.test.default@}


The formula interface converts data into the appropriate contingency matrix for use
with the default method.

@O ../R/multiCA.R @{
#'@@rdname multiCA.test
#'@@method multiCA.test formula
#'@@param formula a formula of the form \code{outcome ~ group} where \code{outcome} is a factor representing the cateogrical outcome and \code{group} is the grouping variableover which the trend is tested.
#'@@param data  an optional matrix or data frame containing the variables in the formula \code{formula}. By default the variables are taken from \code{environment(formula).}
#'@@param subset	an optional vector specifying a subset of observations to be used.
#'@@param na.action	a function which indicates what should happen when the data contain NAs. Defaults to getOption("na.action").
#'@@param weights an integer-valued variable representing the number of times each \code{outcome} - \code{group} combination was observed.
#'@@export

multiCA.test.formula <- function(formula, data, subset, na.action,  weights, ...){
    if (missing(formula) || (length(formula) != 3L) || (length(attr(terms(formula[-2L]), 
        "term.labels")) != 1L)) 
        stop("'formula' missing or incorrect")
    m <- match.call(expand.dots = FALSE)
    if (is.matrix(eval(m$data, parent.frame()))) 
        m$data <- as.data.frame(data)
    m[[1L]] <- quote(stats::model.frame)
    m$... <- NULL
    mf <- eval(m, parent.frame())
    responsevar <- attr(attr(mf, "terms"), "response")
    response <- mf[[responsevar]]
    weightvar <- which(names(mf)=="(weights)")
    w <- if(length(weightvar) > 0)  mf[[weightvar]] else rep(1L, nrow(mf))
    g <- factor(mf[,-c(responsevar, weightvar)])

    tab <- xtabs(w ~ response + g)
    multiCA.test(tab, ...)
}
@| multiCA.test.formula@}

\section{Files}

@f

\section{Macros}

@m

\section{Identifiers}

@u

\end{document}
