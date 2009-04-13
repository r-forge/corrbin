\name{mc.est}
\alias{mc.est}
\title{Distribution of the number of responses assuming marginal compatibility.}
\description{
The \code{mc.est} function estimates the distribution of the number of
responses in a cluster under the assumption of marginal compatibility: information
from all cluster sizes is pooled. The estimation is performed independently for
each treatment group.
}
\usage{
mc.est(cbdata)
}
\arguments{
  \item{cbdata}{a \code{\link{CBData}} object}
}
\details{
The EM algorithm given by Stefanescu and Turnbull (2003) is used.
}
\value{
  A data frame giving the estimated pdf for each treatment and clustersize. It
  has the following columns:
  \item{Prob}{numeric, the probability of \code{NResp} responses in a cluster
    of size \code{ClusterSize} in group \code{Trt}}.  The probabilities add up to 1 for each 
    \code{Trt}/\code{ClusterSize} combination.
  \item{Trt}{factor, the treatment group}
  \item{ClusterSize}{numeric, the cluster size}
  \item{NResp}{numeric, the number of responses}

}
\references{Stefanescu, C. & Turnbull, B. W. (2003) Likelihood inference for exchangeable binary data with varying cluster sizes. 
\emph{Biometrics}, 59, 18-24 }
\author{Aniko Szabo}
\examples{
data(shelltox)
sh.mc <- mc.est(shelltox)

library(lattice)
xyplot(Prob~NResp|factor(ClusterSize), groups=Trt, data=sh.mc, subset=ClusterSize>0, 
    type="l", as.table=TRUE, auto.key=list(columns=4, lines=TRUE, points=FALSE),
    xlab="Number of responses", ylab="Probability P(R=r|N=n)")
}
\keyword{ nonparametric }
\keyword{ models }
