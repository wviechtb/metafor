vcov.matreg <- function(object, ...) {

   mstyle <- .get.mstyle()

   .chkclass(class(object), must="matreg")

   out <- object$vb
   return(out)

}
