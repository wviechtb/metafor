.onAttach <- function(libname, pkgname) {
   ver <- "2.5-20"
   loadmsg <- paste0("Loading 'metafor' package (version ", ver, "). For an overview \nand introduction to the package please type: help(metafor).")
   packageStartupMessage(loadmsg, domain=NULL, appendLF=TRUE)
}
