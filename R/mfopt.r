setmfopt <- function(...) {

   mstyle <- .get.mstyle()

   mfopts <- getOption("metafor")

   if (is.null(mfopts) || !is.list(mfopts)) {
      options("metafor" = list(space=TRUE))
      mfopts <- getOption("metafor")
   }

   newopts <- list(...)

   for (opt in names(newopts)) {
      if (opt == "space" && !is.null(newopts[[opt]]) && !is.logical(newopts[[opt]]))
         stop(mstyle$stop("'space' must be a logical."))
      if (opt == "digits" && !is.null(newopts[[opt]]) && !is.vector(newopts[[opt]], mode="numeric"))
         stop(mstyle$stop("'digits' must be a numeric vector."))
      if (opt == "style" && !is.logical(newopts[[opt]]) && !is.null(newopts[[opt]]) && !is.list(newopts[[opt]]))
         stop(mstyle$stop("'style' must be a list."))
      if (opt == "theme" && !is.null(newopts[[opt]]) && !is.element(newopts[[opt]], c("default", "light", "dark", "auto", "custom", "default2", "light2", "dark2", "auto2", "custom2")))
         stop(mstyle$stop("'theme' must be either 'default(2)', 'light(2)', 'dark(2)', 'auto(2)', or 'custom(2)'."))
      if (opt == "fg" && !is.null(newopts[[opt]]) && !is.character(newopts[[opt]]))
         stop(mstyle$stop("'fg' must be a character string."))
      if (opt == "bg" && !is.null(newopts[[opt]]) && !is.character(newopts[[opt]]))
         stop(mstyle$stop("'bg' must be a character string."))
      mfopts[[opt]] <- newopts[[opt]]
   }

   options("metafor" = mfopts)

}

getmfopt <- function(x, default=NULL) {

   opt <- getOption("metafor")

   if (!missing(x)) {
      x <- as.character(substitute(x))
      opt <- opt[[x]]
   }

   if (is.null(opt)) {
      return(default)
   } else {
      return(opt)
   }

}
