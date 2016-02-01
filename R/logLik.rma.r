logLik.rma <- function(object, REML, ...) {

   if (!inherits(object, "rma"))
      stop("Argument 'object' must be an object of class \"rma\".")

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
   attr(val, "nobs") <- object$k.eff - ifelse(REML, 1, 0) * object$p.eff
   attr(val, "df")   <- object$parms

   class(val) <- "logLik"
   return(val)

}
