.First.lib <- function(libname, pkgname){
  library.dynam("ReprodCalcs", pkgname)
}

.Last.lib <- function(){
  library.dynam.unload("ReprodCalcs")
}
  