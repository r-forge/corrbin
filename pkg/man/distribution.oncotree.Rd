\name{distribution.oncotree}
\alias{distribution.oncotree}
\alias{marginal.distr}
\title{Find the event distribution defined by an oncogenetic tree}
\description{
 \code{distribution.oncotree} calculates the joint distribution
 of the events defined by the tree, while \code{marginal.distr}
 calculates the marginal probability of occurrence of each event.
}
\usage{
   distribution.oncotree(otree, with.probs = TRUE, with.errors=FALSE,
          edge.weights=if (with.errors) "estimated" else "observed")
   marginal.distr(otree, with.errors = TRUE,
          edge.weights=if (with.errors) "estimated" else "observed")
}

\arguments{
  \item{otree}{An object of class \code{oncotree}. }
  \item{with.probs}{A logical value specifying if only the set of
  possible outcomes should be returned (if TRUE), or the associated
  probabilities of occurrence as well.}
  \item{with.errors}{A logical value specifying whether false
  positive and negative error rates should be incorporated into
  the distribution.} 
  \item{edge.weights}{A choice of whether the observed or estimated
  edge transition probabilities should be used in the calculation
  of probabilities. See \code{\link{oncotree.fit}} for explanation
  of the difference. By default, estimated edge transition probabilies
 if \code{with.errors=TRUE} and the observed ones if 
 \code{with.errors=FALSE}.}
}
\value{
  For \code{distribution.oncotree}: a data frame each row of which
  gives a possible outcome. 
 % \describe{
 %   \item{Columns 1--\code{otree$nmut}}{0-1 values for the presence/absence
 %   of each event jointly defining the outcome}
 %   \item{Prob}{the probability of the outcome; present only if 
 %   \code{withprobs=TRUE}.}
 % }
  
  For \code{marginal.distr}: a named numeric vector - the names
  are the event names (+ \sQuote{Root}) and the values are the
  corresponding marginal probability of occurrence.   
}
\author{Aniko Szabo}
\seealso{\code{\link{oncotree.fit}}}
\examples{
   data(ov.cgh)
   ov.tree <- oncotree.fit(ov.cgh)
   
   #joint distribution
   jj <- distribution.oncotree(ov.tree, edge.weights="obs")
   head(jj)
   jj.eps <- distribution.oncotree(ov.tree, with.errors=TRUE)
   head(jj.eps)
  
   #marginal distribution
   marginal.distr(ov.tree, with.error=FALSE)
   #marginal distribution calculated from the joint
   apply(jj[1:ov.tree$nmut], 2, function(x){sum(x*jj$Prob)})
   
   ##Same with errors incorporated
   #marginal distribution
   marginal.distr(ov.tree, with.error=TRUE)
   #marginal distribution calculated from the joint
   apply(jj.eps[1:ov.tree$nmut], 2, function(x){sum(x*jj.eps$Prob)})
   
}
\keyword{models}
\keyword{distribution}
