confint.rma.glmm <- function(object, parm, level, digits, transf, targs, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   if (!inherits(object, "rma.glmm"))
      stop(mstyle$stop("Argument 'object' must be an object of class \"rma.glmm\"."))

   stop(mstyle$stop("Method not yet implemented for objects of class \"rma.glmm\". Sorry!"))

}
