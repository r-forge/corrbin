
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


unwrap.CMData <- function(cmdata){
  #unwrap Freq variable
  freqs <- rep(1:nrow(cmdata), cmdata$Freq)
  cm1 <- cmdata[freqs,]
  cm1$Freq <- NULL
  
  #create ID
  cm1$ID <- factor(1:nrow(cm1))
  
  ncat <- attr(cmdata, "ncat")
  nrespvars <- paste("NResp", 1:ncat, sep=".")
  
  #reshape to have one row per category within cluster
  cm2 <- reshape(cm1, direction="long", varying=list(nrespvars), v.names="Count",
                 idvar="ID", timevar="Outcome", times=1:ncat)
  
  #unwrap categories
  counts <- rep(1:nrow(cm2), cm2$Count)
  res <- cm2[counts,]
  res$Count <- NULL
  class(res) <- "data.frame"
  res <- res[order(res$ID),c("Trt","ID","ClusterSize","Outcome")]
  rownames(res) <- NULL

  res
  }
