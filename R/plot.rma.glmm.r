plot.rma.glmm <- function(x, qqplot=FALSE, ...) {

   #########################################################################

   mstyle <- .get.mstyle("crayon" %in% .packages())

   .chkclass(class(x), must="rma.glmm", notav="rma.glmm")

}
