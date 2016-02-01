# Note: Works with "robust.rma" objects.

coef.rma <- function(object, ...) {

   if (!inherits(object, "rma"))
      stop("Argument 'object' must be an object of class \"rma\".")

   coefs <- c(object$b)
   names(coefs) <- rownames(object$b)
   return(coefs)

}
