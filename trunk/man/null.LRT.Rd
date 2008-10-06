\name{null.LRT}
\alias{null.LRT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
null.LRT(cbdata, control = list())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cbdata}{ ~~Describe \code{cbdata} here~~ }
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
function (cbdata, control = list()) 
{
    a <- with(cbdata, aggregate(Freq, list(ClusterSize = ClusterSize, 
        NResp = NResp), sum))
    names(a)[names(a) == "x"] <- "Freq"
    a$ClusterSize <- as.integer(as.character(a$ClusterSize))
    a$NResp <- as.integer(as.character(a$NResp))
    a$Trt <- 1
    b <- reprod.est(a)
    b <- merge(cbdata, b, all.x = T, by = c("ClusterSize", "NResp"))
    ll0 <- with(b, sum(Freq * log(Prob)))
    res <- mix.reprod.mle(cbdata, control = control)
    ll1 <- res$loglik
    lrt <- 2 * (ll1 - ll0)
    attr(lrt, "ll0") <- ll0
    attr(lrt, "ll1") <- ll1
    lrt
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
