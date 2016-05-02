
#'@keywords internal

.multiCA.test <- function(x, scores, outcomes){
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

  #individual tests
  CAT <- X[nonz]^2 / (pdot[nonz] * (1-pdot[nonz])) / s2 
  CAT.p.value <- pchisq(CAT, df=1, lower.tail=FALSE)
  
  #overall test
  if (full || sum(pdot) >= 1){
    Tt <- ( sum(X[nonz]^2 / pdot[nonz])) / s2
  } else {
    Tt <- (sum(X)^2 / (1-sum(pdot)) + sum(X[nonz]^2 / pdot[nonz])) / s2
  }

  df <- length(outcomes) - full
  p.value <- pchisq(Tt, df=df, lower.tail=FALSE)

  res <- list(statistic = Tt, parameter = df, p.value = p.value, 
              indiv.statistics = CAT, indiv.p.value = CAT.p.value)
  return(res)
}

#'@rdname multiCA.test
#'@method multiCA.test default
#'@param scores numeric vector of the same length as the number of ordered groups. Defaults to linearly increasing values
#'@param outcomes integer or character vector defining the set of outcomes (by row index or row name) over which the trend should be tested. Defaults to all outcomes.
#'@param p.adjust.method character string defining the correction method for individual outcome p-values. Defaults to "closed.set" when \code{length(outcomes)<=3}, and "holm-schaffer" otherwise.
#'@export

multiCA.test.default <- function(x, scores=1:ncol(x), outcomes=1:nrow(x),
  p.adjust.method=c("none","closed.set","holm-schaffer"),...){
  if (!is.matrix(x)) {
    cat(str(x))
    stop("x should be a two-dimensional matrix")
}
  if (length(scores) != ncol(x)) stop("The length of the score vector should equal the number of columns of x")

  testres <- .multiCA.test(x=x, scores=scores, outcomes=outcomes)
 
  Tt <- c(W = testres$statistic)
  df <- c(df = testres$parameter)

  p.value <- testres$p.value
  null.value <- 0
  names(null.value) <- sprintf("slope for outcomes %s", deparse(substitute(outcomes)))

  res <- list(statistic = Tt, parameter = df, p.value = p.value, 
              method="Multinomial Cochran-Armitage trend test",
              alternative="two.sided",
              null.value=null.value,
              data.name = deparse(substitute(x)))
  class(res) <- "htest"

  
    if (missing(p.adjust.method)){
      if (length(outcomes)<=3) p.adjust.method <- "closed.set"
      else p.adjust.method <- "holm-schaffer"
    } else {
      p.adjust.method <- match.arg(p.adjust.method)
    }

    full.set <- (length(outcomes) == nrow(x)) 
    if (p.adjust.method=="none") {
      indiv.res <- testres$indiv.p.value
    } else if (p.adjust.method=="closed.set") {
      
        mytest <- function(hypotheses){
          .multiCA.test(x, scores, hypotheses)$p.value
        }
        indiv.res <- .p.adjust.closed(mytest, outcomes, remove=full.set)  
      
    } else if (p.adjust.method=="holm-schaffer") {
      
          s <- seq_along(testres$indiv.p.value)
          if (full.set) s[2] <- 3
          o <- order(testres$indiv.p.value)
          ro <- order(o)
          indiv.res <- pmin(1, cummax((length(outcomes) - s + 1L) * testres$indiv.p.value[o]))[ro]
      
    } 
  

  return(list(overall = res, individual = indiv.res))  
}

#'@rdname multiCA.test
#'@method multiCA.test formula
#'@param formula a formula of the form \code{outcome ~ group} where \code{outcome} is a factor representing the cateogrical outcome and \code{group} is the grouping variableover which the trend is tested.
#'@param data  an optional matrix or data frame containing the variables in the formula \code{formula}. By default the variables are taken from \code{environment(formula).}
#'@param subset an optional vector specifying a subset of observations to be used.
#'@param na.action      a function which indicates what should happen when the data contain NAs. Defaults to getOption("na.action").
#'@param weights an integer-valued variable representing the number of times each \code{outcome} - \code{group} combination was observed.
#'@export

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

#' @importFrom bitops bitAnd
#' @keywords internal
.bit2boolean <- function (x, N) 
{
  base <- 2^(1:N - 1)
  bitAnd(x, base) != 0
}

#' @param test function that performs the local test. The function should accept a subvector of the hypotheses argument as input, and return a p-value.
#' @param  hypotheses identifiers of the collection of elementary hypotheses. 
#' @param remove logical indicator of whether hypotheses of length N-1 should be removed
#' @param  ...  additional parameters to the 'test' function
#' @return  numeric vector of adjusted p-values for each hypothesis
#' @keywords internal
.p.adjust.closed <- function (test, hypotheses, remove=FALSE, ...) 
{
  N <- length(hypotheses)
  Nmax <- log2(.Machine$integer.max + 1)
  if (N > Nmax) 
    stop("no more than ", Nmax, " hypotheses supported in full closed testing.\n Use a shortcut-based test.")
  closure <- 1:(2^N - 1)
  base <- 2^(1:N - 1)
  offspring <- function(x) {
    res <- bitAnd(x, closure)
    res[res != 0]
  }
  lengths <- rowSums(sapply(base, function(bs) bitAnd(closure, bs) != 0))
  
  idx <- sort.list(lengths, decreasing = TRUE)
  closure <- closure[idx]
  lengths <- lengths[idx]
  if (remove)  closure <- closure[lengths != (N-1)]
  
  adjusted <- numeric(2^N - 1)
  for (i in closure) {
    if (adjusted[i] < 1) {
      localtest <- test(hypotheses[.bit2boolean(i,N)], ...)
      if (localtest > adjusted[i]) {
        offs <- offspring(i)
        adjusted[offs] <- pmax(adjusted[offs], localtest)
      }
    }
  }
  
  out <- adjusted[base]
  names(out) <- hypotheses
  return(out)
}

#' Non-centrality parameter for chi-square distribution
#'
#' Calculates the non-centrality parameter for a chi-square distribution that achieves the target power at a given significance level. This is often needed for sample size calculation for chi-square based tests.
#'@param df an integer giving the degrees of freedom of the chi-square variable
#'@param alpha a numeric value giving the significance level of the test
#'@param beta a numeric value giving the desired type II error (1-\code{beta} is the power)
#'@examples
#' cnonct(6, 0.05, 0.2)
#'@export

cnonct <- function(df, alpha, beta){
  crit.value <- qchisq(alpha, df=df, lower.tail=FALSE)
  
  f <- function(ncp){pchisq(crit.value, df=df, ncp=pmax(0,ncp)) - beta}

  res <- uniroot(f, interval=c(0, 100), extendInt="downX")
  res$root
}
