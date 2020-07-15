.onAttach <- function(libname, pkgname) {

   ver <- "2.5-28"
   loadmsg <- paste0("Loading 'metafor' package (version ", ver, "). For an overview \nand introduction to the package please type: help(metafor).")

   options("pboptions" = list(
      type = if (interactive()) "timer" else "none",
      char = "=",
      txt.width = 50,
      gui.width = 300,
      style = 3,
      initial = 0,
      title = "Progress Bar",
      label = "",
      nout = 100L,
      min_time = 2,
      use_lb = FALSE))

   packageStartupMessage(loadmsg, domain=NULL, appendLF=TRUE)

}
