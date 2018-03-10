weights.rma.glmm <- function(object, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   if (!inherits(object, "rma.glmm"))
      stop(mstyle$stop("Argument 'object' must be an object of class \"rma.glmm\"."))

   stop(mstyle$stop("Method not available for objects of class \"rma.glmm\"."))

}
