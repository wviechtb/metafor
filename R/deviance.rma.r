deviance.rma <- function(object, REML, ...) {

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
      val <- object$fit.stats["dev","REML"]
   } else {
      val <- object$fit.stats["dev","ML"]
   }

   return(val)

}
