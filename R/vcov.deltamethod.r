vcov.deltamethod <- function(object, ...) {

   mstyle <- .get.mstyle()

   .chkclass(class(object), must="deltamethod")

   out <- object$vcov
   return(out)

}
