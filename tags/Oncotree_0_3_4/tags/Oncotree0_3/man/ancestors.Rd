\name{ancestors}
\alias{ancestors}
\alias{least.common.ancestor}
\title{Find ancestors within an oncogenetic tree.}
\description{
\code{ancestors} finds all the ancestors of the given vertex within the tree starting from itself up to the root. \code{least.common.ancestor} finds the common ancestor of two vertices that is closest to them (and farthest from the root).
}
\usage{
   ancestors(otree, vertex)
   least.common.ancestor(otree, v1, v2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{otree}{An object of class \code{oncotree}.}
  \item{vertex, v1, v2}{Character values giving the names of the nodes.}
}
\value{
  For \code{ancestors}: a character vector giving the names of the ancestors of \code{vertex}. The first element is \code{vertex}, and the last one is \dQuote{Root}.
  
   For \code{least.common.ancestor}: a character value with the name of the least common ancestor of \code{v1} and \code{v2}.

}
\seealso{\code{\link{oncotree.fit}}}
\examples{
  data(ov.cgh)
  ov.tree <- oncotree.fit(ov.cgh)
  ancestors(ov.tree, "4q-")
  ancestors(ov.tree, "Xp-")
  least.common.ancestor(ov.tree, "4q-","Xp-")  #"5q-"
}
\keyword{models}
