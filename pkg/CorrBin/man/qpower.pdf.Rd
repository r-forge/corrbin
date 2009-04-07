\name{qpower.pdf}
\alias{qpower.pdf}
\alias{betabin.pdf}
\title{Parametric distributions for correlated binary data}
\description{
\code{qpower.pdf} and \code{betabin.pdf} calculate the probability distribution
function for the number of responses in a cluster of the q-power and beta-binomial 
distributions, respectively.
}
\usage{
qpower.pdf(p, rho, n)
betabin.pdf(p, rho, n)
}
\arguments{
  \item{p}{numeric, the probability of success. }
  \item{rho}{numeric between 0 and 1 inclusive, the within-cluster correlation. }
  \item{n}{integer, cluster size.}
}
\details{
The pdf of the q-power distribution is
 \deqn{P(X=x) = C(n,x)\sum_{k=0}^x (-1)^kC(x,k)q^((n-x+k)^g),}{P(X=x) = \binom{n}{x}\sum_{k=0}^x (-1)^k\binom{x}{k}q^{(n-x+k)^\gamma},}
\eqn{x=0,..,n}{x=0,\ldots,n}, where \eqn{q=1-p}{q=1-p}, and the intra-cluster correlation
 \deqn{rho = (q^(2^g)-q^2)/(q(1-q)).}{\rho = \frac{(q^{2^\gamma}-q^2}{q(1-q)}.}

 The pdf of the beta-binomial distribution is
  \deqn{P(X=x) = C(n,x) B(a+x,n+b-x)/B(a,b),}{P(X=x) = \binom{n}{x} \frac{B(\alpha+x, n+\beta-x)}{B(\alpha,\beta)},}
\eqn{x=0,..,n}{x=0,\ldots,n}, where  \eqn{a=p(1-rho)/rho}{\alpha= p\frac{1-\rho}{\rho}}, and 
\eqn{b=(1-p)(1-rho)/rho}{\alpha= (1-p)\frac{1-\rho}{\rho}}.
}
\value{
a numeric vector of length \eqn{n+1}{n+1} giving the value P(X=x) for \eqn{x=0,..,n}{x=0,\ldots,n}.
}
\references{
Kuk, A. A (2004) litter-based approach to risk assessement in developmental toxicity studies via a power family of completely monotone functions \emph{Applied Statistics}, 52, 51-61.

Williams, D. A. (1975) The Analysis of Binary Responses from Toxicological Experiments Involving Reproduction and Teratogenicity \emph{Biometrics}, 31, 949-952.}
\author{Aniko Szabo, aszabo@mcw.edu }
\seealso{\code{\link{ran.CBData}} for generating an entire dataset using these functions}
\examples{
#the distributions have quite different shapes
#with q-power assigning more weight to the "all affected" event than other distributions
plot(0:10, betabin.pdf(0.3, 0.4, 10), type="o", ylim=c(0,0.34), ylab="Density", xlab="Number of responses out of 10")
lines(0:10, qpower.pdf(0.3, 0.4, 10), type="o", col="red")
}
\keyword{ distribution }
