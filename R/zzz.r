.onAttach <- function(libname, pkgname) {
   loadmsg <- "Loading 'metafor' package (version 2.1-0). For an overview \nand introduction to the package please type: help(metafor)."
   packageStartupMessage(loadmsg, domain=NULL, appendLF=TRUE)
}
