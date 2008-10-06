\name{mix.mc.mle}
\alias{mix.mc.mle}
\title{}
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mix.reprod.mle(cbdata, turn = 1, control = list())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cbdata}{ ~~Describe \code{cbdata} here~~ }
  \item{turn}{ ~~Describe \code{turn} here~~ }
  \item{control}{ ~~Describe \code{control} here~~ }
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\references{ ~put references to the literature/web site here ~ }
\author{ ~~who you are~~ }
\note{ ~~further notes~~ 

 ~Make other sections like Warning with \section{Warning }{....} ~
}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (cbdata, turn = 1, control = list()) 
{
    attach(control)
    on.exit(detach(control))
    tab <- xtabs(Freq ~ factor(ClusterSize, levels = 1:max(ClusterSize)) + 
        factor(NResp, levels = 0:max(ClusterSize)) + Trt, data = cbdata)
    size <- dim(tab)[1]
    ntrt <- dim(tab)[3]
    ntot <- sum(tab)
    storage.mode(tab) <- "double"
    Q <- array(0, dim = rep(size + 1, ntrt))
    storage.mode(Q) <- "double"
    S <- DownUpMatrix(size, ntrt, turn)
    storage.mode(S) <- "integer"
    if ((start == "H0") & (method == "EM")) {
        warning("The EM algorithm can only use 'start=uniform'. Switching options.")
        start <- "uniform"
    }
    if (start == "H0") {
        const.row <- matrix(0:size, nr = size + 1, nc = ntrt)
        Q[const.row + 1] <- 1/(size + 1)
    }
    else {
        Q[S + 1] <- 1/(nrow(S))
    }
    res <- switch(method, EM = .Call("MixReprodQ", Q, S, tab, 
        as.integer(max.iter), as.double(eps), as.integer(verbose)), 
        ISDM = .Call("ReprodISDM", Q, S, tab, as.integer(max.iter), 
            as.integer(max.directions), as.double(eps), as.integer(verbose)))
    names(res) <- c("MLest", "Q", "D", "loglik", "converge")
    dimnames(res$MLest) <- list(NResp = 0:size, ClusterSize = 1:size, 
        Trt = 1:ntrt)
    res$MLest <- as.data.frame.table(res$MLest)
    names(res$MLest) <- c("NResp", "ClusterSize", "Trt", "Prob")
    res$MLest$NResp <- as.numeric(as.character(res$MLest$NResp))
    res$MLest$ClusterSize <- as.numeric(as.character(res$MLest$ClusterSize))
    res$MLest <- subset(res$MLest, NResp <= ClusterSize)
    levels(res$MLest$Trt) <- levels(cbdata$Trt)
    names(res$converge) <- c("rel.error", "n.iter")
    res
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
