df.residual.rma <- function(object, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   if (!inherits(object, "rma"))
      stop(mstyle$stop("Argument 'object' must be an object of class \"rma\"."))

   df.resid <- object$k.eff - object$p.eff

   return(df.resid)

}
