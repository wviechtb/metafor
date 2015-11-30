.onAttach <- function(libname, pkgname) {
   loadmsg <- "Loading 'metafor' package (version 1.9-9). For an overview \nand introduction to the package please type: help(metafor)."
   packageStartupMessage(loadmsg, domain=NULL, appendLF=TRUE)
}
