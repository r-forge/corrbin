
#'Create a `CMdata' object from a data frame.
#'
#'The \code{CMData} function creates an object of class \dfn{CMData} that is
#'used in further analyses. It identifies the variables that define treatment
#'group, clustersize and the number of responses for each outcome type.
#'
#'@export
#'@param x a data frame with one row representing a cluster or potentially a
#'set of clusters of the same size and number of responses for each outcome
#'@param trt the name of the variable that defines treatment group
#'@param nresp either a character vector with the names or a numeric vector with indices
#'of the variables that define the number of responses in
#'the cluster for each outcome type. If \code{clustersize} is \code{NULL}, then it will be
#'calculated assuming that the \code{nresp} vector contains all the possible outcomes.
#'If \code{clustersize} is given, then an additional category is created for the excess cluster members.
#'@param clustersize the name or index of the variable that defines cluster size, or \code{NULL}. If \code{NULL},
#'its value will be calculated by adding the counts from the \code{nresp} variables. If defined,
#'an additional response type will be created for the excess cluster members.
#'@param freq the name or index of the variable that defines the number of clusters
#'represented by the data row. If \code{NULL}, then each row is assumed to
#'correspond to one cluster.
#'@return A data frame with each row representing all the clusters with the
#'same trt/size/number of responses, and standardized variable names:
#'@return \item{Trt}{factor, the treatment group}
#'@return \item{ClusterSize}{numeric, the cluster size}
#'@return \item{NResp.1--NResp.K+1}{numeric, the number of responses for each of the K+1 outcome types}
#'@return \item{Freq}{numeric, number of clusters with the same values}
#'@author Aniko Szabo
#'@seealso \code{\link{read.CMData}} for creating a \code{CMData} object
#'directly from a file.
#'@keywords classes manip
#'@examples
#'
#'data(dehp)
#'dehp <- CMData(dehp, trt="Trt", nresp=c("NResp.1","NResp.2","NResp.3"))
#'str(dehp)
#'

CMData <- function(x, trt, nresp, clustersize=NULL, freq=NULL){
  if (!is.data.frame(x)) stop("x has to be a data frame")
  nms <- names(x)
  K <- if (is.null(clustersize)) length(nresp)-1 else length(nresp)
  
  process.var <- function(var){
    if (is.character(var)){
       if (var %in% nms) res <- x[[var]]
       else stop(paste("Variable '", var, "' not found"))
    }
    else {
      if (is.numeric(var)){
         if (var %in% seq(along=nms)) res <- x[[var]]
         else stop(paste("Column", var, " not found"))
      }
      else stop(paste("Invalid variable specification:",var))
    }
  }
  trtvar <- factor(process.var(trt))
  nrespvar <- sapply(nresp, process.var)
  if (is.null(freq)) freqvar <- rep(1, nrow(x))
  else freqvar <- process.var(freq)
  
  if (!is.null(clustersize)){
     csvar <- process.var(clustersize) # read cluster sizes
     nrespvar <- cbind(nrespvar, csvar - rowSums(nrespvar))  # calculate last category
  }
  else {
    csvar <- rowSums(nrespvar) #calculate sample sizes
  }
  colnames(nrespvar) <- 1:(K+1)
     
  
  d <- data.frame(Trt=trtvar, ClusterSize=csvar, NResp=nrespvar, Freq=freqvar)
  nrespnames <- grep("NResp", names(d), value=TRUE)
  d <- aggregate(d$Freq, d[,c("Trt", "ClusterSize", nrespnames)], sum)
  names(d)[length(names(d))] <- "Freq"
  
  attr(d, "ncat") <- K+1
  class(d) <- c("CMData", "data.frame")
  d}

#'Read data from external file into a CMData object
#'
#'A convenience function to read data from specially structured file directly
#'into a \code{CMData} object. There are two basic data format options:  etiher the counts of responses of all categories are given (and the
#'cluster size is the sum of these counts), or  the total cluster size is given with the counts of all but one category.
#'The first column should always give the treatment group, then either the counts for each category (first option, chosen by setting 
#'\code{with.clustersize = FALSE}), or the size of the cluster followed by the counts for all but one category (second option,
#'chosen by setting \code{with.clustersize = TRUE}). Optionally, a last column could
#'give the number of times the given combination occurs in the data.
#'
#'@export
#'@param file name of file with data. The data in the file should be structured as described above.
#'@param with.clustersize logical indicator of whether a cluster size variable is present
#'in the file
#'@param with.freq logical indicator of whether a frequency variable is present
#'in the file
#'@param ... additonal arguments passed to \code{\link[utils]{read.table}}
#'@return a \code{CMData} object
#'@author Aniko Szabo
#'@seealso \code{\link{CMData}}
#'@keywords IO file
#'

read.CMData <- function(file, with.clustersize=TRUE, with.freq=TRUE, ...){
  d <- read.table(file, ...)
  K <- ncol(d) - with.freq - 2  #subtracting Trt & either ClusterSize or last category column
  
  if (with.clustersize){
    d2 <- CMData(d, trt=1, clustersize=2, nresp=3:(K+2), freq=if (with.freq) "Freq" else NULL) 
  }
  else {
    d2 <- CMData(d, trt=1, nresp=2:(K+2), freq=if (with.freq) "Freq" else NULL)
  }
  d2}

#'\code{unwrap.CMData} is a utility function that reformats a CMData object so
#'that each row is one observation (instead of one or more clusters). A new
#'`ID' variable is added to indicate clusters. This form can be useful for
#'setting up the data for a different package.
#'
#'@rdname unwrap
#'@method unwrap CMData
#'@S3method unwrap CMData
#'@export
#'@return For \code{unwrap.CMData}: a data frame with one row for each cluster element (having a multinomial
#'outcome) with the following standardized column names
#'@return \item{Trt}{factor, the treatment group}
#'@return \item{ClusterSize}{numeric, the cluster size}
#'@return \item{ID}{factor, each level representing a different cluster}
#'@return \item{Resp}{numeric with integer values giving the response type of the cluster
#'element}
#'@examples
#'
#'data(dehp)
#'dehp.long <- unwrap(dehp)
#'head(dehp.long)
#'


unwrap.CMData <- function(object,...){
  #unwrap Freq variable
  freqs <- rep(1:nrow(object), object$Freq)
  cm1 <- object[freqs,]
  cm1$Freq <- NULL
  
  #create ID
  cm1$ID <- factor(1:nrow(cm1))
  
  ncat <- attr(object, "ncat")
  nrespvars <- paste("NResp", 1:ncat, sep=".")
  
  #reshape to have one row per category within cluster
  cm2 <- reshape(cm1, direction="long", varying=list(nrespvars), v.names="Count",
                 idvar="ID", timevar="Resp", times=1:ncat)
  
  #unwrap categories
  counts <- rep(1:nrow(cm2), cm2$Count)
  res <- cm2[counts,]
  res$Count <- NULL
  class(res) <- "data.frame"
  res <- res[order(res$ID),c("Trt","ID","ClusterSize","Resp")]
  rownames(res) <- NULL

  res
  }
