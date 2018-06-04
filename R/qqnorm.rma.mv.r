qqnorm.rma.mv <- function(y, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   if (!inherits(y, "rma.mv"))
      stop(mstyle$stop("Argument 'y' must be an object of class \"rma.mv\"."))

   stop(mstyle$stop("Method not available for objects of class \"rma.mv\"."))

}
