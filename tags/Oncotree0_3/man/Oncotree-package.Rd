\name{Oncotree-package}
\alias{Oncotree-package}
\alias{Oncotree}
\docType{package}
\title{
Constructing and evaluating oncogenetic trees
}
\description{
Oncogenetic trees are directed tree structures that  model the process of
occurrence of genetic alterations during carcinogenesis.
}
\details{
\tabular{ll}{
Package: \tab Oncotree\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2007-01-30\cr
License: \tab GPL 2 or newer\cr
}

A \bold{pure oncogenetic tree} is a directed rooted tree \emph{T} with a
probability \eqn{\pi (e)} attached to each edge \emph{e} such that
for every vertex there is a unique directed path from the root to
it along the edges of the tree. This tree
generates observations on the presence/absence of genetic events the following
way: each edge \emph{e} is independently retained with probability \eqn{\pi
(e)}; the set of vertices that are still reachable from the root 
gives the set of the observed genetic evets.

To describe random deviations from the pure tree model an error model
is added. 

\strong{Error model}
\enumerate{
\item The tumor develops according to the pure oncogenetic tree
model
\item The presence/absence of each alteration is independently
measured
\item If the alteration is present it is not observed with
probability \eqn{\epsilon_-}.

If the alteration is absent it is observed with
probability \eqn{\epsilon_+}.
}
}
\author{
Lisa Pappas, Aniko Szabo

Maintainer: Aniko Szabo <aszabo@mcw.edu>
}
\references{
 [1] Desper R., Jiang F., Kallioniemi O.P., Moch H., Papadimitriou C.H., and
Sch\"affer  A.A. (1999) Inferring tree models for oncogenesis from comparative genome
  hybridization data. \emph{Journal of Computational Biology}. \bold{6}m 37--51.
 [2] Szabo, A. and Boucher, K. (2002) Estimating an oncogenetic
 tree when false negative and positives are present. \emph{Mathematical Biosciences},
 \bold{176}/2, 219--236.
}
\keyword{ package }

\examples{
  data(ov.cgh)
  ov.tree <- oncotree.fit(ov.cgh)
  plot(ov.tree, edge.weights="estimated")
}
