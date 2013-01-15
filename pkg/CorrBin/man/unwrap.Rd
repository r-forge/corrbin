\name{unwrap}
\alias{unwrap}
\alias{unwrap.CBData}
\title{Unwrap a CBData object}
\usage{
  unwrap(object, ...)

  \method{unwrap}{CBData} (object, ...)
}
\arguments{
  \item{object}{a \code{\link{CBData}} object}

  \item{\dots}{other potential arguments; not currently
  used}
}
\value{
  A data frame with one row for each cluster element
  (having a binary outcome) with the following standardized
  column names

  \item{Trt}{factor, the treatment group}

  \item{ClusterSize}{numeric, the cluster size}

  \item{ID}{factor, each level representing a different
  cluster}

  \item{Resp}{numeric with 0/1 values, giving the response
  of the cluster element}
}
\description{
  \code{unwrap.CBData} is a utility function that reformats
  a CBData object so that each row is one observation
  (instead of one or more clusters). A new `ID' variable is
  added to indicate clusters. This form can be useful for
  setting up the data for a different package.
}
\examples{
data(shelltox)
ush <- unwrap(shelltox)
head(ush)
}
\author{
  Aniko Szabo
}
\keyword{manip}

