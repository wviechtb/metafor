# Note: Works with "robust.rma" objects.

df.residual.rma <- function(object, ...) {

   if (!is.element("rma", class(object)))
      stop("Argument 'object' must be an object of class \"rma\".")

   df.resid <- object$k.eff - object$p.eff

   return(df.resid)

}
