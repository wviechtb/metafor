qqnorm.rma.glmm <- function(y, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   if (!inherits(y, "rma.glmm"))
      stop(mstyle$stop("Argument 'y' must be an object of class \"rma.glmm\"."))

   stop(mstyle$stop("Method not yet implemented for objects of class \"rma.glmm\". Sorry!"))

}
