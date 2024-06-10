.onAttach <- function(libname, pkgname) {

   ver <- "4.7-15"

   loadmsg <- paste0("\nLoading the 'metafor' package (version ", ver, "). For an\nintroduction to the package please type: help(metafor)\n")

   installed.ver <- as.numeric(strsplit(gsub("-", ".", ver, fixed=TRUE), ".", fixed=TRUE)[[1]])

   # set default options

   mfopts <- getOption("metafor")

   if (is.null(mfopts) || !is.list(mfopts)) {
      options("metafor" = list(check=TRUE, silent=FALSE, space=TRUE, theme="default"))
   } else {
      if (is.null(mfopts$check))
         mfopts$check <- TRUE
      if (is.null(mfopts$silent))
         mfopts$silent <- FALSE
      if (is.null(mfopts$space))
         mfopts$space <- TRUE
      if (is.null(mfopts$theme))
         mfopts$theme <- "default"
      options("metafor" = mfopts)
   }

   # only run version check in an interactive session and if METAFOR_VERSION_CHECK is not FALSE

   verchk <- tolower(Sys.getenv("METAFOR_VERSION_CHECK")) # "" if unset

   checkopt <- getOption("metafor")$check

   if (!is.null(checkopt)) {
      if (is.logical(checkopt) && isFALSE(checkopt))
         verchk <- "false"
      if (is.character(checkopt) && isTRUE(checkopt == "devel"))
         verchk <- "devel"
   }

   if (interactive() && verchk != "false") {

      #print("Version check ...")

      if (isTRUE(verchk == "devel")) {

         # pull version number from GitHub

         tmp <- suppressWarnings(try(readLines("https://raw.githubusercontent.com/wviechtb/metafor/master/DESCRIPTION", n=2), silent=TRUE))

         if (!inherits(tmp, "try-error") && length(tmp) == 2L) {
            available.ver <- tmp[2]
            if (!is.na(available.ver) && length(available.ver) != 0L)
               available.ver <- substr(available.ver, 10, nchar(available.ver)) # strip 'Version: ' part
         }

      } else {

         # pull version number from CRAN

         tmp <- suppressWarnings(try(readLines("https://cran.r-project.org/web/packages/metafor/index.html"), silent=TRUE))

         if (!inherits(tmp, "try-error")) {
            available.ver <- tmp[grep("Version:", tmp, fixed=TRUE) + 1]
            if (!is.na(available.ver) && length(available.ver) != 0L)
               available.ver <- substr(available.ver, 5, nchar(available.ver)-5) # strip <td> and </td>
         }

      }

      if (!inherits(tmp, "try-error")) {
         save.ver <- available.ver # need this below is message
         available.ver <- as.numeric(strsplit(gsub("-", ".", available.ver), ".", fixed=TRUE)[[1]])
         installed.ver <- 100000 * installed.ver[1] + 1000 * installed.ver[2] + installed.ver[3]
         available.ver <- 100000 * available.ver[1] + 1000 * available.ver[2] + available.ver[3]
         if (isTRUE(installed.ver < available.ver)) {
            loadmsg <- paste0(loadmsg, "\nAn updated version of the package (version ", save.ver, ") is available!\nTo update to this version type: ")
            if (isTRUE(verchk == "devel")) {
               loadmsg <- paste0(loadmsg, "remotes::install_github(\"wviechtb/metafor\")\n")
            } else {
               loadmsg <- paste0(loadmsg, "install.packages(\"metafor\")\n")
            }
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

   if (isFALSE(getOption("metafor")$silent))
      packageStartupMessage(loadmsg, domain=NULL, appendLF=TRUE)

}

.metafor <- new.env()
