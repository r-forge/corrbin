\name{mixControl}
\alias{mixControl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mixControl(method = c("ISDM", "EM"), eps = 0.005, max.iter = 5000, max.directions = 0, start = ifelse(method == "ISDM", "H0", "uniform"), verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{method}{ ~~Describe \code{method} here~~ }
  \item{eps}{ ~~Describe \code{eps} here~~ }
  \item{max.iter}{ ~~Describe \code{max.iter} here~~ }
  \item{max.directions}{ ~~Describe \code{max.directions} here~~ }
  \item{start}{ ~~Describe \code{start} here~~ }
  \item{verbose}{ ~~Describe \code{verbose} here~~ }
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
function (method = c("ISDM", "EM"), eps = 0.005, max.iter = 5000, 
    max.directions = 0, start = ifelse(method == "ISDM", "H0", 
        "uniform"), verbose = FALSE) 
{
    method <- match.arg(method)
    start <- match.arg(start, c("uniform", "H0"))
    list(method = match.arg(method), eps = eps, max.iter = max.iter, 
        max.directions = max.directions, start = start, verbose = verbose)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
