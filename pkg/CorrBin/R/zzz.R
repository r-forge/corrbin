.First.lib <- function(libname, pkgname){
  library.dynam("CorrBin", pkgname, libname)
}

.Last.lib <- function(){
  library.dynam.unload("CorrBin")
}
  