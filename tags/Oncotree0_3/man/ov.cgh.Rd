\name{ov.cgh}
\alias{ov.cgh}
\docType{data}
\title{Ovarian cancer CGH data}
\description{
This is a data set obtained using the comparative genomic hybridization
technique (CGH) on samples from papillary serous cystadenocarcinoma
of the ovary. Only the seven most commonly occurring events are
given.}

\usage{data(ov.cgh)}
\format{
  A data frame with 87 observations on the following 7 variables.
  \describe{
    \item{\code{8q+}}{a 0/1 indicator of the presence of the \sQuote{8q+} event}
    \item{\code{3q+}}{a 0/1 indicator of the presence of the \sQuote{3q+} event}
    \item{\code{5q-}}{a 0/1 indicator of the presence of the \sQuote{5q-} event}
    \item{\code{4q-}}{a 0/1 indicator of the presence of the \sQuote{4q-} event}
    \item{\code{8p-}}{a 0/1 indicator of the presence of the \sQuote{8p-} event}
    \item{\code{1q+}}{a 0/1 indicator of the presence of the \sQuote{1q+} event}
    \item{\code{Xp-}}{a 0/1 indicator of the presence of the \sQuote{Xp-} event}
  }
}
\details{
 The CGH technique uses fluorescent staining to
detect abnormal (increased or decreased) number of DNA copies.
Often the results are reported as a gain or loss on a certain arm,
without further distinction for specific regions. It is common to 
denote a change in DNA copy number on a specific chromosome arm by prefixing
a \dQuote{-} sign for decrease and a \dQuote{+} for increase. Thus, 
say, -3q denotes abnormally low DNA copy number on the q arm of
the 3rd chromosome.
}
\source{
  \url{http://www.ncbi.nlm.nih.gov/sky/}
}
\examples{
  data(ov.cgh)
  heatmap(data.matrix(ov.cgh), Colv=NA, scale="none", col=c("gray90","red"))
}
\keyword{datasets}
