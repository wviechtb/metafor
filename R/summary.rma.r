summary.rma <- function(object, digits, ...) {

   mstyle <- .get.mstyle()

   .chkclass(class(object), must="rma")

   if (missing(digits)) {
      digits <- .get.digits(xdigits=object$digits, dmiss=TRUE)
   } else {
      digits <- .get.digits(digits=digits, xdigits=object$digits, dmiss=FALSE)
   }

   object$digits <- digits

   class(object) <- c("summary.rma", class(object))
   return(object)

}
