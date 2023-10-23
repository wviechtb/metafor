logLik.rma <- function(object, REML, ...) {

   mstyle <- .get.mstyle()

   .chkclass(class(object), must="rma")

   # in case something like logLik(res1, res2) is used

   if (!missing(REML) && inherits(REML, "rma"))
      REML <- NULL

   if (missing(REML) || is.null(REML)) {
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
