.onAttach <- function(libname, pkgname) {

   ver <- "2.5-85"

   loadmsg <- paste0("Loading 'metafor' package (version ", ver, "). For an overview \nand introduction to the package please type: help(metafor)")

   inst.ver <- as.numeric(strsplit(gsub("-", ".", ver, fixed=TRUE), ".", fixed=TRUE)[[1]])

   if (inst.ver[2] %% 2 == 0) {

      # only run version check if a non-devel version is installed

      cran.ver <- suppressWarnings(try(readLines("https://raw.githubusercontent.com/wviechtb/metafor/master/CRAN_version", n=1), silent=TRUE))

      if (!inherits(cran.ver, "try-error")) {
         save.ver <- cran.ver
         cran.ver <- as.numeric(strsplit(gsub("-", ".", cran.ver), ".", fixed=TRUE)[[1]])
         inst.ver <- 100000 * inst.ver[1] + 1000 * inst.ver[2] + inst.ver[3]
         cran.ver <- 100000 * cran.ver[1] + 1000 * cran.ver[2] + cran.ver[3]
         if (inst.ver < cran.ver)
            loadmsg <- paste0(loadmsg, "\n\nAn updated version of the package (version ", save.ver, ") is available!\nTo update to this version type: install.packages(\"metafor\")")
      }

   }

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
