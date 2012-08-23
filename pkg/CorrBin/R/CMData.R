
CMData <- function(x, trt, clustersize, nresp, freq=NULL){
  if (!is.data.frame(x)) stop("x has to be a data frame")
  nms <- names(x)
  K <- length(nresp)
  
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
  csvar <- process.var(clustersize)
  nrespvar <- sapply(nresp, process.var)
  if (K > 1) colnames(nrespvar) <- 1:K
  if (is.null(freq)) freqvar <- rep(1, nrow(x))
  else freqvar <- process.var(freq)
  
  d <- data.frame(Trt=trtvar, ClusterSize=csvar, NResp=nrespvar, Freq=freqvar)
  nrespnames <- grep("NResp", names(d), value=TRUE)
  d <- aggregate(d$Freq, d[,c("Trt", "ClusterSize", nrespnames)], sum)
  names(d)[length(names(d))] <- "Freq"
  
  attr(d, "K") <- K
  class(d) <- c("CMData", "data.frame")
  d}

read.CMData <- function(file, with.clustersize=TRUE, with.freq=TRUE, ...){
  d <- read.table(file, ...)
  K <- ncol(d) - with.freq - 2  #subtracting Trt & either ClusterSize or last category column
  nrespvars <- paste("NResp", 1:K, sep=".")
  if (with.clustersize) 
    names(d) <- c("Trt","ClusterSize", nrespvars, if (with.freq) "Freq"))
  else {
   names(d) <- c("Trt", nrespvars, "LastCat", if (with.freq) "Freq"))
   d$ClusterSize <- rowSums(d[, c(nrespvars, "LastCat")]
   d$LastCat <- NULL
  }
  d <- CMData(d, trt="Trt", clustersize="ClusterSize", nresp=nrespvars, 
              freq=if (with.freq) "Freq" else NULL)
  d}


unwrap.CMData <- function(cmdata){
  #unwrap Freq variable
  freqs <- rep(1:nrow(cmdata), cmdata$Freq)
  cm1 <- cmdata[freqs,]
  cm1$Freq <- NULL
  
  #create ID
  cm1$ID <- factor(1:nrow(cm1))
  
  #create last category count
  nrespvars <- grep("NResp", names(cm1), value=TRUE)
  cm1$LastCat <- cm1$ClusterSize - rowSums(cm1[,nrespvars])
  K <- length(nrespvars)
  
  #reshape to have one row per category within cluster
  cm2 <- reshape(cm1, direction="long", varying=list(c(nrespvars, "LastCat")), v.names="Count",
                 idvar="ID", timevar="Outcome", times=1:(K+1))
  
  #unwrap categories
  counts <- rep(1:nrow(cm2), cm2$Count)
  res <- cm2[counts,]
  cm2$Count <- NULL

  res[order(res$ID),]
  }
