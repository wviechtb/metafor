qqnorm.rma.mv <- function(y, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   .chkclass(class(y), must="rma.mv", notav="rma.mv")

}
