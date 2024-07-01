nobs.rma <- function(object, ...) {

   mstyle <- .get.mstyle()

   .chkclass(class(object), must="rma")

   #n.obs <- object$k.eff - ifelse(object$method == "REML", 1, 0) * object$p.eff
   n.obs <- object$k

   return(n.obs)

}
