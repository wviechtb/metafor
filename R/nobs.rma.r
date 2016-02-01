nobs.rma <- function(object, ...) {

   if (!inherits(object, "rma"))
      stop("Argument 'object' must be an object of class \"rma\".")

   n.obs <- object$k.eff - ifelse(object$method == "REML", 1, 0) * object$p.eff

   return(n.obs)

}
