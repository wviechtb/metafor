plot.rma.glmm <- function(x, qqplot=FALSE, ...) {

   #########################################################################

   mstyle <- .get.mstyle("crayon" %in% .packages())

   if (!inherits(x, "rma.glmm"))
      stop(mstyle$stop("Argument 'x' must be an object of class \"rma.glmm\"."))

   stop(mstyle$stop("Method not available for objects of class \"rma.glmm\"."))

}
