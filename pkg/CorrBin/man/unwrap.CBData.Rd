\name{unwrap.CBData}
\alias{unwrap.CBData}
\title{Unwrap a CBData object}
\description{
\code{unwrap.CBData} is a utility function that reformats a CBData object so that
each row is one observation (instead of one or more clusters). A new `ID' variable is added
to indicate clusters. This form can be useful for setting up the data for a different package.
}
\usage{
unwrap.CBData(cbdata)
}
\arguments{
  \item{cbdata}{a \code{\link{CBData}} object}
}
\value{
A data frame with one row for each cluster element (having a binary outcome) with the
following standardized column names  
  \item{Trt}{factor, the treatment group}
  \item{ClusterSize}{numeric, the cluster size}
  \item{ID}{factor, each level representing a different cluster}
  \item{Resp}{numeric with 0/1 values, giving the response of the cluster element}
 
}
\author{ Aniko Szabo}
\examples{
data(shelltox)
ush <- unwrap.CBData(shelltox)
head(ush)
}
\keyword{ manip }

