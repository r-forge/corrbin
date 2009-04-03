\name{bootstrap}
\alias{bootstrap.oncotree}
\alias{print.boottree}
\alias{plot.boottree}
\title{Bootstrap an oncogenetic tree to assess stability}
\description{\code{bootstrap.oncotree} provides a set of resampling
based estimates of the oncogenetic tree. Both a parametric and 
non-parametric approach is available. The \code{print} and 
\code{plot} methods provide interfaces for printing a summary and
plotting the resulting set of trees.}

\usage{
   bootstrap.oncotree(otree, R, type = c("nonparametric", "parametric"))
   \method{print}{boottree}(x, \dots)
   \method{plot}{boottree}(x, minfreq=NULL, minprop=NULL, nboots=NULL, draw.orig=TRUE, draw.consensus=TRUE, 
	                         fix.nodes=FALSE, ask=(prod(par("mfrow"))<ntrees)&&dev.interactive(), \dots)
}
\arguments{
  \item{otree}{An object of class \code{oncotree}.}
  \item{R}{The number of bootstrap replicates.}
  \item{type}{The type of bootstrap - see Details for explanations.}
  \item{x}{An object of class \code{boottree} - the output of \code{bootstrap.oncotree}}
  \item{minfreq}{A lower limit on the occurrence frequency of the tree in \dQuote{boottree} for plotting. By default, all unique trees are plotted, which can lead to a large number of plots.}
  \item{minprop}{A lower limit on the occurrence proportion of the tree in \dQuote{boottree} for plotting.}
  \item{nboots}{A lower limit on the number of bootstrapped trees plotted.}
  \item{draw.orig}{logical; if TRUE the original tree is plotted.}
  \item{draw.consensus}{logical; if TRUE the consensus tree is plotted (see Details).}
  \item{fix.nodes}{logical; if TRUE, the nodes for all trees are kept in the same position. If \code{node.coords}
	is passed as an argument to \code{\link{plot.oncotree}}, then those coordinates are used for all trees, otherwise
	the coordinates computed for the original tree are used.}
  \item{ask}{logical; if TRUE, the user is \emph{ask}ed before each plot, see \code{\link{par}}(ask=.).}
  \item{\dots}{Ignored for \code{print}. Passed to \code{\link{plot.oncotree}}
  for the \code{plot} method.}
}
\details{
 \bold{Parametric bootstrap:} This approach assumes that the model is correct. Based on \code{otree}, a random data set is generated R times using \code{\link{generate.data}}. An oncogenetic tree is fitted to each of these random data sets. 
 
 \bold{Non-parametric bootstrap:} The samples (rows) from the data associated with the tree are resampled with replacement R times, each time obtaining a data set with the same sample size. An oncogenetic tree is fitted to each of these resampled data sets. 
 
 For both approaches, a \emph{consensus tree} that assigns to each vertex the parent that occurs most frequently in the bootstrapped trees, is also computed.
}
\value{
  For \code{bootstrap.oncotree}: an object of class \code{boottree} with the following components:
  \item{original}{The \code{parent} component of the original tree (\code{otree}).}
  \item{consensus}{A numeric vector with the \code{parent$parent.num} component of the consesus tree - this defines the tree structure uniquely.}
  \item{parent.freq}{A matrix giving the number of trees with each possible child-parent edge. The rows correspond to children while the column to parents.}
	\item{tree.list}{A data frame with each row representing a unique tree obtained during the bootstrap.  The \sQuote{Tree} variable contains the \code{parent$parent.num} component of the tree (each pasted into one dot-separated string), while the \sQuote{Freq} variable gives the frequency of the tree among the R bootstrap replicates.}
  \item{type}{A character value with the type of the bootstrap performed.}
  
  For \code{print.boottree}:
  the original object is returned invisibly. It prints a summary showing the number of replicates, the number of unique trees found, and the number of times that the original tree was obtained.
  
For \code{plot.oncotree}:
  nothing is returned. It is used for its side effect of producing a sequence of plots of the bootstrapped trees. Specifically, it plots the original tree (if \code{draw.orig=TRUE}), the consensus tree (if \code{draw.consensus=TRUE}), and then the other trees by frequency of occurrence. To limit the number of bootstrapped trees plotted, specify exactly one of \code{minfreq}, \code{minprop} or \code{nboots}. By default, if the session is interactive, the user is asked for confirmation before each new tree is drawn. To avoid this, either use \code{ask=FALSE} in the function call, or set up a layout that fits all the trees. 
}
\author{Lisa Pappas, Aniko Szabo}
\seealso{\code{\link{oncotree.fit}}}
\examples{
   data(ov.cgh)
   ov.tree <- oncotree.fit(ov.cgh)
   set.seed(43636)
   ov.b1 <- bootstrap.oncotree(ov.tree, R=100, type="parametric")
   ov.b1
   opar <- par(mfrow=c(3,2), mar=c(2,0,0,0))
   plot(ov.b1, nboots=4)
   plot(ov.b1, nboots=4, fix.nodes=TRUE)
   par(opar)
}
\keyword{models}
\keyword{nonparametric}
