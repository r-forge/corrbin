
.First.lib <- function(libname, pkgname){
  library.dynam("CorrBin", pkgname, libname)
}

.Last.lib <- function(libpath){
  library.dynam.unload("CorrBin",libpath)
}
