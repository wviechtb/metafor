summary.rma <- function(object, digits, showfit=TRUE, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   if (!inherits(object, "rma"))
      stop(mstyle$stop("Argument 'object' must be an object of class \"rma\"."))

   if (missing(digits))
      digits <- object$digits

   object$digits <- digits

   class(object) <- c("summary.rma", class(object))
   return(object)

}
