deviance.rma <- function(object, REML, ...) {

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
      val <- object$fit.stats["dev","REML"]
   } else {
      val <- object$fit.stats["dev","ML"]
   }

   return(val)

}
