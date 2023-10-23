df.residual.rma <- function(object, ...) {

   mstyle <- .get.mstyle()

   .chkclass(class(object), must="rma")

   df.resid <- object$k.eff - object$p.eff

   return(df.resid)

}
