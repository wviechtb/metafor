.onAttach <- function(libname, pkgname) {
   ver <- "2.2-18"
   loadmsg <- paste0("Loading 'metafor' package (version ", ver, "). For an overview \nand introduction to the package please type: help(metafor).")
   packageStartupMessage(loadmsg, domain=NULL, appendLF=TRUE)
}
