############################################################################

as.data.frame.confint.rma <- function(x, ...) {

   .chkclass(class(x), must="confint.rma")

   ddd <- list(...)

   .chkdots(ddd, c("fixed", "random"))

   if (is.null(ddd$fixed)) {
      fixed <- is.element("fixed", names(x))
   } else {
      fixed <- ddd$fixed
   }

   if (is.null(ddd$random)) {
      random <- is.element("random", names(x))
   } else {
      random <- ddd$random
   }

   if (fixed) {
      df <- x$fixed
   } else {
      df <- NULL
   }

   if (random && is.element("random", names(x)))
      df <- rbind(df, x$random)

   return(df)

}

as.data.frame.list.confint.rma <- function(x, ...) {

   .chkclass(class(x), must="list.confint.rma")

   x$digits <- NULL # remove digits elements

   df <- lapply(x, as.data.frame)
   df <- do.call(rbind, df)
   return(df)

}

############################################################################
