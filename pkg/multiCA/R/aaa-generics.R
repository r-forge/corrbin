
#'Multinomial Cochran-Armitage trend test
#'
#'The \code{multiCA.test} performs a multinomial generalization of the 
#' Cochran-Armitage trend test.
#'
#'
#'@export
#'@param x a two-dimensional matrix or a formula
#'@param \dots other arguments 
#'@return an object of class "htest" with the results of the test
#'@author Aniko Szabo
#'@references Szabo, A. (2016) Test for trend with a multinomial outcome.  
#'@keywords nonparametric 
#'@examples
#'
#'data(stroke)
#'## using formula interface
#'multiCA.test(Type ~ Year, weights=Freq, data=stroke)
#'
#'## using matrix interface and testing only the first 3 outcomes
#'strk.mat <- xtabs(Freq ~ Type + Year, data=stroke)
#'multiCA.test(strk.mat, outcomes=1:3)
#'
#'@name multiCA.test

multiCA.test <- function(x,...) UseMethod("multiCA.test")
 
