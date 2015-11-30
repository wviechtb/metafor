# Note: Works with "robust.rma" objects.

summary.rma <- function(object, digits, showfit=TRUE, ...) {

   if (!is.element("rma", class(object)))
      stop("Argument 'object' must be an object of class \"rma\".")

   if (missing(digits))
      digits <- object$digits

   object$digits <- digits

   class(object) <- c("summary.rma", class(object))
   return(object)

}
