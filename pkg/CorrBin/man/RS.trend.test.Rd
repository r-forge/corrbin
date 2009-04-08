\name{RS.trend.test}
\alias{RS.trend.test}
\title{Rao-Scott trend test}
\description{
\code{RS.trend.test} implements the Rao-Scott adjusted Cochran-Armitage test for 
linear increasing trend with correlated data. 
}
\usage{
RS.trend.test(cbdata)
}
\arguments{
  \item{cbdata}{a \code{\link{CBData}} object}
}
\details{
 The test is based on calculating a \dfn{design effect} for each cluster by
 dividing the observed variability by the one expected under independence. The
 number of responses and the cluster size are then divided by the design effect,
 and a Cochran-Armitage type test statistic is computed based on these adjusted
 values.
 
 The implementation aims for testing for \emph{increasing} trend, and a one-sided
 p-value is reported. The test statistic is asymptotically normally distributed,
 and a two-sided p-value can be easily computed if needed.
}
\value{
A list with components
  \item{statistic }{numeric, the value of the test statistic}
  \item{p.val }{numeric, asymptotic one-sided p-value of the test}
}
\references{
Rao, J. N. K. & Scott, A. J. A  (1992) Simple Method for the Analysis of Clustered Data \emph{Biometrics}, 48, 577-586.
 }
\author{Aniko Szabo, aszabo@mcw.edu}
\seealso{ \code{\link{SO.trend.test}}, \code{\link{GEE.trend.test}} for alternative tests;
\code{\link{CBData}} for constructing a CBData object.}
\examples{
data(shelltox)
RS.trend.test(shelltox)
}
\keyword{ htest }
\keyword{ nonparametric }
