.First.lib <- function(libname, pkgname){
  library.dynam("ReprodCalcs", pkgname, libname)
}

.Last.lib <- function(){
  library.dynam.unload("ReprodCalcs")
}
  