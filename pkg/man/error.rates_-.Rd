\name{error.rates<-}
\alias{error.rates<-}
\title{Set the error rates of an oncotree manually}
\description{
 Allows to set the false positive and false negative error rate associated with an object of class \code{oncotree} to values other than those found by the optimization in \code{\link{oncotree.fit}}. The estimated edge transition probabilities are updated appropriately.
}
\usage{
error.rates(x) <- value
}
\arguments{
  \item{x}{An object of class \code{oncotree}.}
  \item{value}{A numeric vector of length 2. The false positive error rate will be set to \code{value[1]}, while the false negative error rate to \code{value[2]}.}
}
\seealso{\code{\link{oncotree.fit}}}
\examples{
  data(ov.cgh)
  ov.tree <- oncotree.fit(ov.cgh)
  ov.tree
  error.rates(ov.tree) <- c(0,0)
  ov.tree
}
\keyword{models}
