confint.rma.glmm <- function(object, parm, level, digits, transf, targs, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   .chkclass(class(object), must="rma.glmm", notav="rma.glmm")

}
