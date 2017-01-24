summary.rma <- function(object, digits, showfit=TRUE, ...) {

   if (!inherits(object, "rma"))
      stop("Argument 'object' must be an object of class \"rma\".")

   if (missing(digits))
      digits <- object$digits

   object$digits <- digits

   class(object) <- c("summary.rma", class(object))
   return(object)

}
