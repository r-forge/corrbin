\name{unwrap}
\alias{unwrap}
\alias{unwrap.CBData}
\alias{unwrap.CMData}
\title{Unwrap a clustered object}
\usage{
unwrap(object, ...)

\method{unwrap}{CBData}(object, ...)

\method{unwrap}{CMData}(object, ...)
}
\arguments{
  \item{object}{a \code{\link{CBData}} object}

  \item{\dots}{other potential arguments; not currently
  used}
}
\value{
For \code{uwrap.CBData}: a data frame with one row for each
cluster element (having a binary outcome) with the
following standardized column names

\item{Trt}{factor, the treatment group}

\item{ClusterSize}{numeric, the cluster size}

\item{ID}{factor, each level representing a different
cluster}

\item{Resp}{numeric with 0/1 values, giving the response of
the cluster element}

For \code{unwrap.CMData}: a data frame with one row for
each cluster element (having a multinomial outcome) with
the following standardized column names

\item{Trt}{factor, the treatment group}

\item{ClusterSize}{numeric, the cluster size}

\item{ID}{factor, each level representing a different
cluster}

\item{Resp}{numeric with integer values giving the response
type of the cluster element}
}
\description{
\code{unwrap.CBData} is a utility function that reformats a
CBData object so that each row is one observation (instead
of one or more clusters). A new `ID' variable is added to
indicate clusters. This form can be useful for setting up
the data for a different package.

\code{unwrap.CMData} is a utility function that reformats a
CMData object so that each row is one observation (instead
of one or more clusters). A new `ID' variable is added to
indicate clusters. This form can be useful for setting up
the data for a different package.
}
\examples{
data(shelltox)
ush <- unwrap(shelltox)
head(ush)
data(dehp)
dehp.long <- unwrap(dehp)
head(dehp.long)
}
\author{
Aniko Szabo
}
\keyword{manip}

