.onAttach <- function(libname, pkgname) {

   ver <- "3.1-34"

   loadmsg <- paste0("\nLoading the 'metafor' package (version ", ver, "). For an\nintroduction to the package please type: help(metafor)\n")

   inst.ver <- as.numeric(strsplit(gsub("-", ".", ver, fixed=TRUE), ".", fixed=TRUE)[[1]])

   # only run version check in an interactive session and if METAFOR_VERSION_CHECK is not FALSE

   verchk <- tolower(Sys.getenv("METAFOR_VERSION_CHECK"))

   if (interactive() && verchk != "false") {

      # pull version number from CRAN page

      tmp <- suppressWarnings(try(readLines("https://cran.r-project.org/web/packages/metafor/index.html"), silent=TRUE))

      # or pull version number from github

      # tmp <- suppressWarnings(try(readLines("https://raw.githubusercontent.com/wviechtb/metafor/master/CRAN_version", n=1), silent=TRUE))

      if (!inherits(tmp, "try-error")) {
         cran.ver <- tmp[grep("Version:", tmp, fixed=TRUE) + 1]
         if (!is.na(cran.ver) && length(cran.ver) != 0L) {
            cran.ver <- substr(cran.ver, 5, nchar(cran.ver)-5) # strip <td> and </td>
            save.ver <- cran.ver # need this below is message
            cran.ver <- as.numeric(strsplit(gsub("-", ".", cran.ver), ".", fixed=TRUE)[[1]])
            inst.ver <- 100000 * inst.ver[1] + 1000 * inst.ver[2] + inst.ver[3]
            cran.ver <- 100000 * cran.ver[1] + 1000 * cran.ver[2] + cran.ver[3]
            if (isTRUE(inst.ver < cran.ver))
               loadmsg <- paste0(loadmsg, "\nAn updated version of the package (version ", save.ver, ") is available!\nTo update to this version type: install.packages(\"metafor\")\n")
         }
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
