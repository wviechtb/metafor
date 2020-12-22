deviance.rma <- function(object, REML, ...) {

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
      val <- object$fit.stats["dev","REML"]
   } else {
      val <- object$fit.stats["dev","ML"]
   }

   return(val)

}
