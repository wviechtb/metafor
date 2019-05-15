.onAttach <- function(libname, pkgname) {
   loadmsg <- "Loading 'metafor' package (version 2.2-0). For an overview \nand introduction to the package please type: help(metafor)."
   packageStartupMessage(loadmsg, domain=NULL, appendLF=TRUE)
}
