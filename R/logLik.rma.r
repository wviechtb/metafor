logLik.rma <- function(object, REML, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   .chkclass(class(object), must="rma")

   if (missing(REML)) {
      if (object$method == "REML") {
         REML <- TRUE
      } else {
         REML <- FALSE
      }
   }

   if (REML) {
      val <- object$fit.stats["ll","REML"]
   } else {
      val <- object$fit.stats["ll","ML"]
   }

   attr(val, "nall") <- object$k.eff
   attr(val, "nobs") <- object$k.eff - ifelse(REML, object$p.eff, 0)
   attr(val, "df")   <- object$parms

   class(val) <- "logLik"
   return(val)

}
