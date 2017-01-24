df.residual.rma <- function(object, ...) {

   if (!inherits(object, "rma"))
      stop("Argument 'object' must be an object of class \"rma\".")

   df.resid <- object$k.eff - object$p.eff

   return(df.resid)

}
