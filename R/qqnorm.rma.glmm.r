qqnorm.rma.glmm <- function(y, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   .chkclass(class(y), must="rma.glmm", notav="rma.glmm")

}
