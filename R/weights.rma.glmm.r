weights.rma.glmm <- function(object, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   .chkclass(class(object), must="rma.glmm", notav="rma.glmm")

}
