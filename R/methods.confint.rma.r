############################################################################

as.data.frame.confint.rma <- function(x, ...) {

   .chkclass(class(x), must="confint.rma")

   ddd <- list(...)

   .chkdots(ddd, c("fixed", "random"))

   fixed  <- .chkddd(ddd$fixed,  is.element("fixed",  names(x)))
   random <- .chkddd(ddd$random, is.element("random", names(x)))

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
