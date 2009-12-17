\name{oncotree}
\alias{oncotree.fit}
\alias{print.oncotree}
\alias{plot.oncotree}
\alias{pstree.oncotree}
\title{Build and display an oncogenetic tree}
\description{
 Build a directed tree structure to model the process of occurrence
 of genetic alterations (events) in carcinogenesis. The model is
 described in more detail in \code{\link{Oncotree-package}}. Methods for 
 printing a short summary, displaying the tree on an R plot, and producing
 latex code for drawing the tree (using the \sQuote{pstricks} and \sQuote{pst-tree} LaTeX packages) 
 are provided.
}
\usage{
  oncotree.fit(dataset, error.fun = function(x, y){sum((x - y)^2)})
  \method{print}{oncotree}(x, \dots)
  \method{plot}{oncotree}(x, edge.weights = c("none", "observed", "estimated"), 
            edge.digits=2, node.coords=NULL, plot=TRUE, cex = par("cex"), 
						col.edge=par("col"), col.text=par("col"), col.weight=par("col"),\dots)
  pstree.oncotree(x, edge.weights=c("none","observed","estimated"), edge.digits=2,
                  shape=c("none","oval", "circle", "triangle", "diamond"),
                  pstree.options=list(arrows="->", treefit="loose", arrowscale="1.5 0.8", nodesep="3pt"))
}


\arguments{
  \item{dataset}{A data frame or a matrix with variable names as 
  a listing of genetic events taking on binary values indicating 
  missing (0) or present (1). Each row is an independent sample.}
  \item{error.fun}{A function of two variables that measures the 
  deviation of the observed marginal frequencies of the events 
  (which will be the first argument in the call) from the estimated ones. 
  The false positive and negative error rates are obtained by  minimizing 
  \code{error.fun}. If \code{error.fun=NULL} is used, the error rates are not estimated.}
  \item{x}{An object of class \code{oncotree}.}
  \item{edge.weights}{Choice of edge weights to show on the plot.}
  \item{edge.digits}{The number of significant digits to use when displaying edge weights.}
  \item{node.coords}{A matrix with node-coordinates or NULL if the coordinates should be
	computed automatically (default).}
  \item{plot}{Logical; indicates whether the tree should be plotted.}
  \item{cex}{Scaling factor for the text in the nodes.}
  \item{col.edge}{color of the tree edges.}
  \item{col.text}{color of the node label.}
  \item{col.weight}{color of the edge weights.}
  \item{\dots}{Ignored for \code{print}. For \code{plot} these can
  be graphical parameters passed to \code{lines} when the edges are
  drawn}
	\item{shape}{The shape of the node in the pst-tree representation.}
	\item{pstree.options}{Additional options for pst-tree. See the pstricks
	documentation for possible values.}
}
\value{
For \code{oncotree.fit}:
  an object of class \code{oncotree} which has components
\item{data}{data frame used, after dropping events with zero observed frequency,
            and adding a column for the artificial \sQuote{Root} node}
\item{nmut}{number of tree nodes: the number genetic events present in data +1 for the \sQuote{Root} node}
\item{parent}{a list containing information about the tree structure with the following components
  \itemize{ 
  \item{child}{a character vector of the event names starting with \sQuote{Root}}
  \item{parent}{a character vector of the names of the parents of \code{child}}
  \item{parent.num}{a numeric vector with column indices corresponding to \code{parent}}
  \item{obs.weights}{raw edge transition probabilites P(child|parent)}
  \item{est.weights}{edge transition probabilities adjusted for the error rates \code{eps}}}}
\item{level}{a numeric vector of the depth of each node in the tree (1 for the root, 2 for its children, etc.)}
\item{numchild}{a numeric vector giving the number of children for each node}
\item{levelnodes}{a numeric vector of the number of nodes found at each level of the tree}
\item{levelgrp}{a character matrix with its rows giving the ordered nodes at each level}
\item{eps}{a numeric vector of length two showing the estimated false positive and negative error rates (if \code{error.fun} is not NULL). Do not modify directly, but rather through \code{\link{error.rates<-}}.}

For \code{print.oncotree}:

  the original object is returned invisibly. It prints a summary showing the number of nodes, the parent-child relationships, and the false positive and negative error rates.
  
For \code{plot.oncotree}:

  a matrix with node-coordinates is returned invisibly. The column names of the matrix are the
	names of the nodes/events (including 'Root'), the rows gives the x- and y-coordinates, respectively.
  This matrix provides a valid input for \code{node.coords}. If \code{plot=TRUE}, a plot of the tree is produced.
  
For \code{pstree.oncotree}:

  a character string with the LaTeX code needed to draw a tree. \code{\\usepackage{pstricks,pst-tree}}
  is required in the preamble of the LaTeX file, and it should be processed through
  a PostScript intermediary (DVIPS or similar) and not through PDFLaTeX.

}

\details{
\sQuote{pst-tree} is a very flexible package, and very detailed formatting of the tree
is possible. \code{pstree.oncotree} provides some default settings for drawing
trees, but they can be easily overridden: most options can be set in 
\code{pstree.options}, while the appearance of the tree nodes can be controlled
by defining a one-parameter \code{\\lab} command that gives the desired appearance. 
For example, if red, non-mathematical test is desired in an oval, you could use
\code{\\newcommand{\\lab}[1]{\\Toval[name=#1]{{\red #1}}}}.
}

\references{Szabo, A. and Boucher, K. (2002) Estimating an oncogenetic tree when false negative and positives are present. Mathematical Biosciences, 176/2, 219-236.}
\author{Lisa Pappas}

\seealso{\code{\link{bootstrap.oncotree}},\code{\link{error.rates<-}},
  \code{\link{generate.data}},\code{\link{ancestors}},\code{\link{distribution.oncotree}}}

\examples{
  data(ov.cgh)
  ov.tree <- oncotree.fit(ov.cgh, error.fun=function(x,y){max(abs(x-y))})
  ov.tree
  nodes <- plot(ov.tree, edge.weights="est")
  #move the Root node to the left
  nodes["x","Root"] <- nodes["x","8q+"]
  plot(ov.tree, node.coords=nodes)
  #output for pstricks+pst-tree
  pstree.oncotree(ov.tree, edge.weights="obs", shape="oval")
}
\keyword{models}
\keyword{tree}
\keyword{graphs}
\keyword{hplot}
