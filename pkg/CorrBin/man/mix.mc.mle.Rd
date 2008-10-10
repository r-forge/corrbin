\name{mix.mc.mle}
\alias{mix.mc.mle}
\title{Order-restricted MLE assuming marginal compatibility}
\description{
  \code{mix.mc.mle} computes the nonparametric maximum likelihood estimate
  of the distribution of the number of responses in a cluster \eqn{P(R=r|n)}{P(R=r|n)}
  under a stochastic ordering constraint. Umbrella ordering can be specified using
  the \code{turn} parameter.}
\usage{
mix.reprod.mle(cbdata, turn = 1, control = mixControl())
}
\arguments{
  \item{cbdata}{an object of class \code{\link{CBData}}.}
  \item{turn}{integer specifying the peak of the umbrella ordering (see Details). The
  default corresponds to a non-decreasing order.}
  \item{control}{a list of control parameters, usually a call to \code{\link{mixControl}}.}
}
\details{
  
}
\value{
  A list with components:
  \item{MLest }{data frame with the maximum likelihood estimates of  P(R_i=r|n)}
  \item{Q }{numeric matrix; estimated weights for the mixing distribution}
  \item{D }{numeric matrix; directional derivative of the log-likelihood}
  \item{loglik }{the achieved value of the log-likelihood}
  \item{converge }{a 2-element vector with the achived relative error and the
  performed number of iterations}
Components \code{Q} and \code{D} are unlikely to be needed for the user. 
}
\references{Szabo A, George EO. (2008) On the Use of Stochastic Ordering to Test for Trend with Clustered Binary Data.}
\author{Aniko Szabo, aszabo@mcw.edu}
\seealso{\code{\link{mixControl}}}
\examples{
  data(shelltox)
  ml <- mix.mc.mle(shelltox, control=mixControl(eps=0.01, method="ISDM")))
  ml$converge
  
  require(lattice)
  xyplot(Prob~NResp|factor(ClusterSize), groups=Trt, data=ml$MLest, type="s",
       panel=panel.superpose, panel.groups=panel.cumsum,
       as.table=T, auto.key=list(columns=4, lines=T, points=F),
			 xlab="Number of responses", ylab="Cumulative Probability R(R>=r|N=n)",
			 ylim=c(0,1.1), main="Stochastically ordered estimates\n with marginal compatibility")
}
\keyword{ nonparametric }
\keyword{ models }
