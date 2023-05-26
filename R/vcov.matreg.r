vcov.matreg <- function(object, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   .chkclass(class(object), must="matreg")

   out <- object$vb
   return(out)

}
