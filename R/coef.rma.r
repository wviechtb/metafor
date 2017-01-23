# Note: Works with "robust.rma" objects.

coef.rma <- function(object, ...) {

   if (!inherits(object, "rma"))
      stop("Argument 'object' must be an object of class \"rma\".")

   coefs <- c(object$beta)
   names(coefs) <- rownames(object$beta)

   if (inherits(object, "rma.ls")) {
      coefs <- list(beta=coefs)
      coefs$alpha <- c(object$alpha)
      names(coefs$alpha) <- rownames(object$alpha)
   }

   return(coefs)

}
