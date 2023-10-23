summary.matreg <- function(object, digits, ...) {

   mstyle <- .get.mstyle()

   .chkclass(class(object), must="matreg")

   if (missing(digits)) {
      digits <- .get.digits(xdigits=object$digits, dmiss=TRUE)
   } else {
      digits <- .get.digits(digits=digits, xdigits=object$digits, dmiss=FALSE)
   }

   object$digits <- digits

   class(object) <- c("summary.matreg", class(object))
   return(object)

}
